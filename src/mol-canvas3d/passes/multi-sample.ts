/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { PostprocessingPass } from './postprocessing';
import { DrawPass } from './draw';

import quad_vert from '../../mol-gl/shader/quad.vert';
import compose_frag from '../../mol-gl/shader/compose.frag';

const ComposeSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),
    uWeight: UniformSpec('f'),
};

type ComposeRenderable = ComputeRenderable<Values<typeof ComposeSchema>>

function getComposeRenderable(ctx: WebGLContext, colorTexture: Texture): ComposeRenderable {
    const values: Values<typeof ComposeSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),
        uWeight: ValueCell.create(1.0),
    };

    const schema = { ...ComposeSchema };
    const shaderCode = ShaderCode('compose', quad_vert, compose_frag);
    const renderItem = createComputeRenderItem(ctx, 'triangles', shaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

export const MultiSampleParams = {
    mode: PD.Select('off', [['off', 'Off'], ['on', 'On'], ['temporal', 'Temporal']]),
    sampleLevel: PD.Numeric(2, { min: 0, max: 5, step: 1 }),
};
export type MultiSampleProps = PD.Values<typeof MultiSampleParams>

export class MultiSamplePass {
    props: MultiSampleProps
    colorTarget: RenderTarget

    private composeTarget: RenderTarget
    private holdTarget: RenderTarget
    private compose: ComposeRenderable

    private sampleIndex = -1
    private currentTime = 0
    private lastRenderTime = 0

    constructor(private webgl: WebGLContext, private camera: Camera, private drawPass: DrawPass, private postprocessing: PostprocessingPass, props: Partial<MultiSampleProps>) {
        const { gl } = webgl;
        this.colorTarget = webgl.createRenderTarget(gl.drawingBufferWidth, gl.drawingBufferHeight);
        this.composeTarget = webgl.createRenderTarget(gl.drawingBufferWidth, gl.drawingBufferHeight);
        this.holdTarget = webgl.createRenderTarget(gl.drawingBufferWidth, gl.drawingBufferHeight);
        this.compose = getComposeRenderable(webgl, drawPass.colorTarget.texture);
        this.props = { ...PD.getDefaultValues(MultiSampleParams), ...props };
    }

    get enabled() {
        if (this.props.mode === 'temporal') {
            if (this.currentTime - this.lastRenderTime > 200) {
                return this.sampleIndex !== -1;
            } else {
                this.sampleIndex = 0;
                return false;
            }
        } else if (this.props.mode === 'on') {
            return true;
        } else {
            return false;
        }
    }

    update(changed: boolean, currentTime: number) {
        if (changed) this.lastRenderTime = currentTime;
        this.currentTime = currentTime;
    }

    setSize(width: number, height: number) {
        this.colorTarget.setSize(width, height);
        this.composeTarget.setSize(width, height);
        this.holdTarget.setSize(width, height);
        ValueCell.update(this.compose.values.uTexSize, Vec2.set(this.compose.values.uTexSize.ref.value, width, height));
    }

    setProps(props: Partial<MultiSampleProps>) {
        if (props.mode !== undefined) this.props.mode = props.mode;
        if (props.sampleLevel !== undefined) this.props.sampleLevel = props.sampleLevel;
    }

    render(toDrawingBuffer: boolean, transparentBackground: boolean) {
        if (this.props.mode === 'temporal') {
            this.renderTemporalMultiSample(toDrawingBuffer, transparentBackground);
        } else {
            this.renderMultiSample(toDrawingBuffer, transparentBackground);
        }
    }

    private renderMultiSample(toDrawingBuffer: boolean, transparentBackground: boolean) {
        const { camera, compose, composeTarget, drawPass, postprocessing, webgl } = this;
        const { gl, state } = webgl;

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = JitterVectors[ Math.max(0, Math.min(this.props.sampleLevel, 5)) ];

        const baseSampleWeight = 1.0 / offsetList.length;
        const roundingRange = 1 / 32;

        camera.viewOffset.enabled = true;
        ValueCell.update(compose.values.tColor, postprocessing.enabled ? postprocessing.target.texture : drawPass.colorTarget.texture);
        compose.update();

        const width = drawPass.colorTarget.getWidth();
        const height = drawPass.colorTarget.getHeight();

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

            // render scene and optionally postprocess
            drawPass.render(false, transparentBackground);
            if (postprocessing.enabled) postprocessing.render(false);

            // compose rendered scene with compose target
            composeTarget.bind();
            gl.viewport(0, 0, width, height);
            state.enable(gl.BLEND);
            state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
            state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
            state.disable(gl.DEPTH_TEST);
            state.disable(gl.SCISSOR_TEST);
            state.depthMask(false);
            if (i === 0) {
                state.clearColor(0, 0, 0, 0);
                gl.clear(gl.COLOR_BUFFER_BIT);
            }
            compose.render();
        }

        ValueCell.update(compose.values.uWeight, 1.0);
        ValueCell.update(compose.values.tColor, composeTarget.texture);
        compose.update();

        if (toDrawingBuffer) {
            webgl.unbindFramebuffer();
        } else {
            this.colorTarget.bind();
        }
        gl.viewport(0, 0, width, height);
        state.disable(gl.BLEND);
        compose.render();

        camera.viewOffset.enabled = false;
        camera.update();
    }

    private renderTemporalMultiSample(toDrawingBuffer: boolean, transparentBackground: boolean) {
        const { camera, compose, composeTarget, holdTarget, postprocessing, drawPass, webgl } = this;
        const { gl, state } = webgl;

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = JitterVectors[ Math.max(0, Math.min(this.props.sampleLevel, 5)) ];

        if (this.sampleIndex === -1) return;
        if (this.sampleIndex >= offsetList.length) {
            this.sampleIndex = -1;
            return;
        }

        const i = this.sampleIndex;

        if (i === 0) {
            drawPass.render(false, transparentBackground);
            if (postprocessing.enabled) postprocessing.render(false);
            ValueCell.update(compose.values.uWeight, 1.0);
            ValueCell.update(compose.values.tColor, postprocessing.enabled ? postprocessing.target.texture : drawPass.colorTarget.texture);
            compose.update();

            holdTarget.bind();
            state.disable(gl.BLEND);
            compose.render();
        }

        const sampleWeight = 1.0 / offsetList.length;

        camera.viewOffset.enabled = true;
        ValueCell.update(compose.values.tColor, postprocessing.enabled ? postprocessing.target.texture : drawPass.colorTarget.texture);
        ValueCell.update(compose.values.uWeight, sampleWeight);
        compose.update();

        const width = drawPass.colorTarget.getWidth();
        const height = drawPass.colorTarget.getHeight();

        // render the scene multiple times, each slightly jitter offset
        // from the last and accumulate the results.
        const numSamplesPerFrame = Math.pow(2, this.props.sampleLevel);
        for (let i = 0; i < numSamplesPerFrame; ++i) {
            const offset = offsetList[this.sampleIndex];
            Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height);
            camera.update();

            // render scene and optionally postprocess
            drawPass.render(false, transparentBackground);
            if (postprocessing.enabled) postprocessing.render(false);

            // compose rendered scene with compose target
            composeTarget.bind();
            state.enable(gl.BLEND);
            state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
            state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
            state.disable(gl.DEPTH_TEST);
            state.disable(gl.SCISSOR_TEST);
            state.depthMask(false);
            if (this.sampleIndex === 0) {
                state.clearColor(0, 0, 0, 0);
                gl.clear(gl.COLOR_BUFFER_BIT);
            }
            compose.render();

            this.sampleIndex += 1;
            if (this.sampleIndex >= offsetList.length ) break;
        }

        const accumulationWeight = this.sampleIndex * sampleWeight;
        if (accumulationWeight > 0) {
            ValueCell.update(compose.values.uWeight, 1.0);
            ValueCell.update(compose.values.tColor, composeTarget.texture);
            compose.update();
            if (toDrawingBuffer) {
                webgl.unbindFramebuffer();
            } else {
                this.colorTarget.bind();
            }
            gl.viewport(0, 0, width, height);
            state.disable(gl.BLEND);
            compose.render();
        }
        if (accumulationWeight < 1.0) {
            ValueCell.update(compose.values.uWeight, 1.0 - accumulationWeight);
            ValueCell.update(compose.values.tColor, holdTarget.texture);
            compose.update();
            if (toDrawingBuffer) {
                webgl.unbindFramebuffer();
            } else {
                this.colorTarget.bind();
            }
            gl.viewport(0, 0, width, height);
            if (accumulationWeight === 0) state.disable(gl.BLEND);
            else state.enable(gl.BLEND);
            compose.render();
        }

        camera.viewOffset.enabled = false;
        camera.update();
        if (this.sampleIndex >= offsetList.length) this.sampleIndex = -1;
    }
}

const JitterVectors = [
    [
        [ 0, 0 ]
    ],
    [
        [ 4, 4 ], [ -4, -4 ]
    ],
    [
        [ -2, -6 ], [ 6, -2 ], [ -6, 2 ], [ 2, 6 ]
    ],
    [
        [ 1, -3 ], [ -1, 3 ], [ 5, 1 ], [ -3, -5 ],
        [ -5, 5 ], [ -7, -1 ], [ 3, 7 ], [ 7, -7 ]
    ],
    [
        [ 1, 1 ], [ -1, -3 ], [ -3, 2 ], [ 4, -1 ],
        [ -5, -2 ], [ 2, 5 ], [ 5, 3 ], [ 3, -5 ],
        [ -2, 6 ], [ 0, -7 ], [ -4, -6 ], [ -6, 4 ],
        [ -8, 0 ], [ 7, -4 ], [ 6, 7 ], [ -7, -8 ]
    ],
    [
        [ -4, -7 ], [ -7, -5 ], [ -3, -5 ], [ -5, -4 ],
        [ -1, -4 ], [ -2, -2 ], [ -6, -1 ], [ -4, 0 ],
        [ -7, 1 ], [ -1, 2 ], [ -6, 3 ], [ -3, 3 ],
        [ -7, 6 ], [ -3, 6 ], [ -5, 7 ], [ -1, 7 ],
        [ 5, -7 ], [ 1, -6 ], [ 6, -5 ], [ 4, -4 ],
        [ 2, -3 ], [ 7, -2 ], [ 1, -1 ], [ 4, -1 ],
        [ 2, 1 ], [ 6, 2 ], [ 0, 4 ], [ 4, 4 ],
        [ 2, 5 ], [ 7, 5 ], [ 5, 6 ], [ 3, 7 ]
    ]
];

JitterVectors.forEach(offsetList => {
    offsetList.forEach(offset => {
        // 0.0625 = 1 / 16
        offset[0] *= 0.0625;
        offset[1] *= 0.0625;
    });
});