/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { RenderContext } from '../util';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { compose_frag } from '../../mol-gl/shader/compose.frag';
import { MarkingProps } from './marking';
import { isTimingMode } from '../../mol-util/debug';
import { clearJitter, getJitterOffsets, getTemporalSamplesPerFrame, setJitter } from './jitter';
import { StereoCamera } from '../camera/stereo';

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
    mode: PD.Select('temporal', [['off', 'Off'], ['on', 'On'], ['temporal', 'Temporal']]),
    sampleLevel: PD.Numeric(2, { min: 0, max: 5, step: 1 }, { description: 'Take level^2 samples.' }),
    reduceFlicker: PD.Boolean(true, { description: 'Reduce flicker in "temporal" mode when marking is disabled (i.e. drawn as part of the scene).' }),
    reuseOcclusion: PD.Boolean(true, { description: 'Reuse occlusion data. It is faster but has some artefacts.' }),
};
export type MultiSampleProps = PD.Values<typeof MultiSampleParams>

type Props = {
    multiSample: MultiSampleProps
    postprocessing: PostprocessingProps
    marking: MarkingProps
    transparentBackground: boolean;
    dpoitIterations: number;
}

/** A stereo view can't be marked from a single camera, there marking stays part of the sample rendering. */
function isSeparateMarking(camera: Camera | StereoCamera): camera is Camera {
    return !StereoCamera.is(camera);
}

export class MultiSamplePass {
    static isEnabled(props: MultiSampleProps) {
        return props.mode !== 'off';
    }

    /** accumulated samples, marking is only blended into it when not rendering to the drawing buffer */
    readonly colorTarget: RenderTarget;

    private composeTarget: RenderTarget;
    private holdTarget: RenderTarget;
    private compose: ComposeRenderable;

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

    getByteCount() {
        return this.colorTarget.getByteCount() + this.composeTarget.getByteCount() + this.holdTarget.getByteCount();
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

    render(sampleIndex: number, ctx: RenderContext, props: Props, toDrawingBuffer: boolean, forceOn: boolean) {
        if (props.multiSample.mode === 'temporal' && !forceOn) {
            return this.renderTemporalMultiSample(sampleIndex, ctx, props, toDrawingBuffer);
        } else {
            this.renderMultiSample(ctx, props, toDrawingBuffer);
            return -2;
        }
    }

    private renderMultiSample(ctx: RenderContext, props: Props, toDrawingBuffer: boolean) {
        const { camera } = ctx;
        const { compose, composeTarget, drawPass, webgl } = this;
        const { gl, state } = webgl;
        if (isTimingMode) webgl.timer.mark('MultiSamplePass.renderMultiSample');

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = getJitterOffsets(props.multiSample.sampleLevel);

        const { x, y, width, height } = camera.viewport;
        const baseSampleWeight = 1.0 / offsetList.length;
        const roundingRange = 1 / 32;

        ValueCell.update(compose.values.tColor, drawPass.getColorTarget(props.postprocessing).texture);
        compose.update();

        // render the scene multiple times, each slightly jitter offset
        // from the last and accumulate the results.
        for (let i = 0; i < offsetList.length; ++i) {
            const offset = offsetList[i];
            setJitter(camera, offset);

            // the theory is that equal weights for each sample lead to an accumulation of rounding
            // errors. The following equation varies the sampleWeight per sample so that it is uniformly
            // distributed across a range of values whose rounding errors cancel each other out.
            const uniformCenteredDistribution = -0.5 + (i + 0.5) / offsetList.length;
            const sampleWeight = baseSampleWeight + roundingRange * uniformCenteredDistribution;
            ValueCell.update(compose.values.uWeight, sampleWeight);

            // render scene
            if (i === 0 || !props.multiSample.reuseOcclusion) {
                drawPass.postprocessing.setOcclusionOffset(0, 0);
            } else {
                drawPass.postprocessing.setOcclusionOffset(
                    offset[0] / width,
                    offset[1] / height
                );
            }
            drawPass.render(ctx, props, false, isSeparateMarking(camera));

            // compose rendered scene with compose target
            composeTarget.bind();
            state.enable(gl.BLEND);
            state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
            state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
            state.disable(gl.DEPTH_TEST);
            state.depthMask(false);
            state.viewport(x, y, width, height);
            state.scissor(x, y, width, height);
            if (i === 0) {
                state.clearColor(0, 0, 0, 0);
                gl.clear(gl.COLOR_BUFFER_BIT);
            }
            compose.render();
        }

        drawPass.postprocessing.setOcclusionOffset(0, 0);

        ValueCell.update(compose.values.uWeight, 1.0);
        ValueCell.update(compose.values.tColor, composeTarget.texture);
        compose.update();

        this.colorTarget.bind();
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        state.disable(gl.BLEND);
        compose.render();

        clearJitter(camera);

        this.drawPass.marking.present(ctx, props, { base: this.colorTarget, toDrawingBuffer, offsets: offsetList, restart: true, samples: offsetList.length, shading: drawPass.getMarkingShading(props.postprocessing) });
        if (isTimingMode) webgl.timer.markEnd('MultiSamplePass.renderMultiSample');
    }

    private renderTemporalMultiSample(sampleIndex: number, ctx: RenderContext, props: Props, toDrawingBuffer: boolean) {
        const { camera } = ctx;
        const { compose, composeTarget, holdTarget, drawPass, webgl } = this;
        const { gl, state } = webgl;
        if (isTimingMode) webgl.timer.mark('MultiSamplePass.renderTemporalMultiSample');

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = getJitterOffsets(props.multiSample.sampleLevel);

        if (sampleIndex === -2 || sampleIndex >= offsetList.length) return -2;

        const { x, y, width, height } = camera.viewport;
        const sampleWeight = 1.0 / offsetList.length;
        const samplesPerFrame = getTemporalSamplesPerFrame(props.multiSample.sampleLevel);
        // the unjittered first frame is also the one rendered while the camera moves
        const firstFrame = sampleIndex === -1;

        if (firstFrame) {
            drawPass.render(ctx, props, false, isSeparateMarking(camera));
            ValueCell.update(compose.values.uWeight, 1.0);
            ValueCell.update(compose.values.tColor, drawPass.getColorTarget(props.postprocessing).texture);
            compose.update();

            holdTarget.bind();
            state.disable(gl.BLEND);
            state.disable(gl.DEPTH_TEST);
            state.depthMask(false);
            state.viewport(x, y, width, height);
            state.scissor(x, y, width, height);
            compose.render();
            sampleIndex += 1;
        } else {
            ValueCell.update(compose.values.tColor, drawPass.getColorTarget(props.postprocessing).texture);
            ValueCell.update(compose.values.uWeight, sampleWeight);
            compose.update();

            // render the scene multiple times, each slightly jitter offset
            // from the last and accumulate the results.
            for (let i = 0; i < samplesPerFrame; ++i) {
                const offset = offsetList[sampleIndex];
                setJitter(camera, offset);

                // render scene
                if (sampleIndex === 0 || !props.multiSample.reuseOcclusion) {
                    drawPass.postprocessing.setOcclusionOffset(0, 0);
                } else {
                    drawPass.postprocessing.setOcclusionOffset(
                        offset[0] / width,
                        offset[1] / height
                    );
                }
                drawPass.render(ctx, props, false, isSeparateMarking(camera));

                // compose rendered scene with compose target
                composeTarget.bind();
                state.enable(gl.BLEND);
                state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
                state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
                state.disable(gl.DEPTH_TEST);
                state.depthMask(false);
                state.viewport(x, y, width, height);
                state.scissor(x, y, width, height);
                if (sampleIndex === 0) {
                    state.clearColor(0, 0, 0, 0);
                    gl.clear(gl.COLOR_BUFFER_BIT);
                }
                compose.render();

                sampleIndex += 1;
                if (sampleIndex >= offsetList.length) break;
            }
        }

        drawPass.postprocessing.setOcclusionOffset(0, 0);

        this.colorTarget.bind();
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

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

        clearJitter(camera);

        this.drawPass.marking.present(ctx, props, { base: this.colorTarget, toDrawingBuffer, offsets: offsetList, restart: firstFrame, samples: firstFrame ? 1 : samplesPerFrame, shading: drawPass.getMarkingShading(props.postprocessing) });
        if (isTimingMode) webgl.timer.markEnd('MultiSamplePass.renderTemporalMultiSample');

        return sampleIndex >= offsetList.length ? -2 : sampleIndex;
    }
}

export class MultiSampleHelper {
    private sampleIndex = -2;

    update(changed: boolean, props: MultiSampleProps) {
        if (changed) this.sampleIndex = -1;
        return props.mode === 'temporal' ? this.sampleIndex !== -2 : false;
    }

    /** Return `true` while more samples are needed */
    render(ctx: RenderContext, props: Props, toDrawingBuffer: boolean, forceOn?: boolean) {
        this.sampleIndex = this.multiSamplePass.render(this.sampleIndex, ctx, props, toDrawingBuffer, !!forceOn);
        return this.sampleIndex < 0;
    }

    constructor(private multiSamplePass: MultiSamplePass) {

    }
}