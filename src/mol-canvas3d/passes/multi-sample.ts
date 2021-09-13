/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { Texture } from '../../mol-gl/webgl/texture';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ValueCell } from '../../mol-util';
import { Vec2 } from '../../mol-math/linear-algebra';
import { ShaderCode } from '../../mol-gl/shader-code';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { createComputeRenderable, ComputeRenderable } from '../../mol-gl/renderable';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Camera } from '../../mol-canvas3d/camera';
import { PostprocessingProps } from './postprocessing';
import { DrawPass } from './draw';
import { Renderer } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { Helper } from '../helper/helper';
import { StereoCamera } from '../camera/stereo';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { compose_frag } from '../../mol-gl/shader/compose.frag';
import { MarkingProps } from './marking';

const ComposeSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),
    uWeight: UniformSpec('f'),
};
const ComposeShaderCode = ShaderCode('compose', quad_vert, compose_frag);
type ComposeRenderable = ComputeRenderable<Values<typeof ComposeSchema>>

function getComposeRenderable(ctx: WebGLContext, colorTexture: Texture): ComposeRenderable {
    const values: Values<typeof ComposeSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),
        uWeight: ValueCell.create(1.0),
    };

    const schema = { ...ComposeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', ComposeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export const MultiSampleParams = {
    mode: PD.Select('off', [['off', 'Off'], ['on', 'On'], ['temporal', 'Temporal']]),
    sampleLevel: PD.Numeric(2, { min: 0, max: 5, step: 1 }),
};
export type MultiSampleProps = PD.Values<typeof MultiSampleParams>

type Props = {
    multiSample: MultiSampleProps
    postprocessing: PostprocessingProps
    marking: MarkingProps
}

export class MultiSamplePass {
    static isEnabled(props: MultiSampleProps) {
        return props.mode !== 'off';
    }

    colorTarget: RenderTarget

    private composeTarget: RenderTarget
    private holdTarget: RenderTarget
    private compose: ComposeRenderable

    constructor(private webgl: WebGLContext, private drawPass: DrawPass) {
        const { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } = webgl.extensions;
        const width = drawPass.colorTarget.getWidth();
        const height = drawPass.colorTarget.getHeight();
        this.colorTarget = webgl.createRenderTarget(width, height, false);
        const type = colorBufferHalfFloat && textureHalfFloat ? 'fp16' :
            colorBufferFloat && textureFloat ? 'float32' : 'uint8';
        this.composeTarget = webgl.createRenderTarget(width, height, false, type);
        this.holdTarget = webgl.createRenderTarget(width, height, false);
        this.compose = getComposeRenderable(webgl, drawPass.colorTarget.texture);
    }

    syncSize() {
        const width = this.drawPass.colorTarget.getWidth();
        const height = this.drawPass.colorTarget.getHeight();

        const [w, h] = this.compose.values.uTexSize.ref.value;
        if (width !== w || height !== h) {
            this.colorTarget.setSize(width, height);
            this.composeTarget.setSize(width, height);
            this.holdTarget.setSize(width, height);
            ValueCell.update(this.compose.values.uTexSize, Vec2.set(this.compose.values.uTexSize.ref.value, width, height));
        }
    }

    render(sampleIndex: number, renderer: Renderer, camera: Camera | StereoCamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean, props: Props) {
        if (props.multiSample.mode === 'temporal') {
            return this.renderTemporalMultiSample(sampleIndex, renderer, camera, scene, helper, toDrawingBuffer, transparentBackground, props);
        } else {
            this.renderMultiSample(renderer, camera, scene, helper, toDrawingBuffer, transparentBackground, props);
            return sampleIndex;
        }
    }

    private bindOutputTarget(toDrawingBuffer: boolean) {
        if (toDrawingBuffer) {
            this.webgl.unbindFramebuffer();
        } else {
            this.colorTarget.bind();
        }
    }

    private renderMultiSample(renderer: Renderer, camera: Camera | StereoCamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean, props: Props) {
        const { compose, composeTarget, drawPass, webgl } = this;
        const { gl, state } = webgl;

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = JitterVectors[Math.max(0, Math.min(props.multiSample.sampleLevel, 5))];

        const { x, y, width, height } = camera.viewport;
        const baseSampleWeight = 1.0 / offsetList.length;
        const roundingRange = 1 / 32;

        camera.viewOffset.enabled = true;
        ValueCell.update(compose.values.tColor, drawPass.getColorTarget(props.postprocessing).texture);
        compose.update();

        // render the scene multiple times, each slightly jitter offset
        // from the last and accumulate the results.
        for (let i = 0; i < offsetList.length; ++i) {
            const offset = offsetList[i];
            Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height);
            camera.update();

            // the theory is that equal weights for each sample lead to an accumulation of rounding
            // errors. The following equation varies the sampleWeight per sample so that it is uniformly
            // distributed across a range of values whose rounding errors cancel each other out.
            const uniformCenteredDistribution = -0.5 + (i + 0.5) / offsetList.length;
            const sampleWeight = baseSampleWeight + roundingRange * uniformCenteredDistribution;
            ValueCell.update(compose.values.uWeight, sampleWeight);

            // render scene
            drawPass.render(renderer, camera, scene, helper, false, transparentBackground, props.postprocessing, props.marking);

            // compose rendered scene with compose target
            composeTarget.bind();
            state.enable(gl.BLEND);
            state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
            state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
            state.disable(gl.DEPTH_TEST);
            state.depthMask(false);
            gl.viewport(x, y, width, height);
            gl.scissor(x, y, width, height);
            if (i === 0) {
                state.clearColor(0, 0, 0, 0);
                gl.clear(gl.COLOR_BUFFER_BIT);
            }
            compose.render();
        }

        ValueCell.update(compose.values.uWeight, 1.0);
        ValueCell.update(compose.values.tColor, composeTarget.texture);
        compose.update();

        this.bindOutputTarget(toDrawingBuffer);
        gl.viewport(x, y, width, height);
        gl.scissor(x, y, width, height);

        state.disable(gl.BLEND);
        compose.render();

        camera.viewOffset.enabled = false;
        camera.update();
    }

    private renderTemporalMultiSample(sampleIndex: number, renderer: Renderer, camera: Camera | StereoCamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean, props: Props) {
        const { compose, composeTarget, holdTarget, drawPass, webgl } = this;
        const { gl, state } = webgl;

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = JitterVectors[Math.max(0, Math.min(props.multiSample.sampleLevel, 5))];

        if (sampleIndex === -2 || sampleIndex >= offsetList.length) return -2;

        const { x, y, width, height } = camera.viewport;
        const sampleWeight = 1.0 / offsetList.length;

        if (sampleIndex === -1) {
            drawPass.render(renderer, camera, scene, helper, false, transparentBackground, props.postprocessing, props.marking);
            ValueCell.update(compose.values.uWeight, 1.0);
            ValueCell.update(compose.values.tColor, drawPass.getColorTarget(props.postprocessing).texture);
            compose.update();

            holdTarget.bind();
            state.disable(gl.BLEND);
            state.disable(gl.DEPTH_TEST);
            state.depthMask(false);
            gl.viewport(x, y, width, height);
            gl.scissor(x, y, width, height);
            compose.render();
            sampleIndex += 1;
        } else {
            camera.viewOffset.enabled = true;
            ValueCell.update(compose.values.tColor, drawPass.getColorTarget(props.postprocessing).texture);
            ValueCell.update(compose.values.uWeight, sampleWeight);
            compose.update();

            // render the scene multiple times, each slightly jitter offset
            // from the last and accumulate the results.
            const numSamplesPerFrame = Math.pow(2, Math.max(0, props.multiSample.sampleLevel - 2));
            for (let i = 0; i < numSamplesPerFrame; ++i) {
                const offset = offsetList[sampleIndex];
                Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height);
                camera.update();

                // render scene
                drawPass.render(renderer, camera, scene, helper, false, transparentBackground, props.postprocessing, props.marking);

                // compose rendered scene with compose target
                composeTarget.bind();
                state.enable(gl.BLEND);
                state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
                state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
                state.disable(gl.DEPTH_TEST);
                state.depthMask(false);
                gl.viewport(x, y, width, height);
                gl.scissor(x, y, width, height);
                if (sampleIndex === 0) {
                    state.clearColor(0, 0, 0, 0);
                    gl.clear(gl.COLOR_BUFFER_BIT);
                }
                compose.render();

                sampleIndex += 1;
                if (sampleIndex >= offsetList.length) break;
            }
        }

        this.bindOutputTarget(toDrawingBuffer);
        gl.viewport(x, y, width, height);
        gl.scissor(x, y, width, height);

        const accumulationWeight = sampleIndex * sampleWeight;
        if (accumulationWeight > 0) {
            ValueCell.update(compose.values.uWeight, 1.0);
            ValueCell.update(compose.values.tColor, composeTarget.texture);
            compose.update();
            state.disable(gl.BLEND);
            compose.render();
        }
        if (accumulationWeight < 1.0) {
            ValueCell.update(compose.values.uWeight, 1.0 - accumulationWeight);
            ValueCell.update(compose.values.tColor, holdTarget.texture);
            compose.update();
            if (accumulationWeight === 0) state.disable(gl.BLEND);
            else state.enable(gl.BLEND);
            compose.render();
        }

        camera.viewOffset.enabled = false;
        camera.update();

        return sampleIndex >= offsetList.length ? -2 : sampleIndex;
    }
}

const JitterVectors = [
    [
        [0, 0]
    ],
    [
        [4, 4], [-4, -4]
    ],
    [
        [-2, -6], [6, -2], [-6, 2], [2, 6]
    ],
    [
        [1, -3], [-1, 3], [5, 1], [-3, -5],
        [-5, 5], [-7, -1], [3, 7], [7, -7]
    ],
    [
        [1, 1], [-1, -3], [-3, 2], [4, -1],
        [-5, -2], [2, 5], [5, 3], [3, -5],
        [-2, 6], [0, -7], [-4, -6], [-6, 4],
        [-8, 0], [7, -4], [6, 7], [-7, -8]
    ],
    [
        [-4, -7], [-7, -5], [-3, -5], [-5, -4],
        [-1, -4], [-2, -2], [-6, -1], [-4, 0],
        [-7, 1], [-1, 2], [-6, 3], [-3, 3],
        [-7, 6], [-3, 6], [-5, 7], [-1, 7],
        [5, -7], [1, -6], [6, -5], [4, -4],
        [2, -3], [7, -2], [1, -1], [4, -1],
        [2, 1], [6, 2], [0, 4], [4, 4],
        [2, 5], [7, 5], [5, 6], [3, 7]
    ]
];

JitterVectors.forEach(offsetList => {
    offsetList.forEach(offset => {
        // 0.0625 = 1 / 16
        offset[0] *= 0.0625;
        offset[1] *= 0.0625;
    });
});

export class MultiSampleHelper {
    private sampleIndex = -2

    update(changed: boolean, props: MultiSampleProps) {
        if (changed) this.sampleIndex = -1;
        return props.mode === 'temporal' ? this.sampleIndex !== -2 : false;
    }

    render(renderer: Renderer, camera: Camera | StereoCamera, scene: Scene, helper: Helper, toDrawingBuffer: boolean, transparentBackground: boolean, props: Props) {
        this.sampleIndex = this.multiSamplePass.render(this.sampleIndex, renderer, camera, scene, helper, toDrawingBuffer, transparentBackground, props);
    }

    constructor(private multiSamplePass: MultiSamplePass) {

    }
}