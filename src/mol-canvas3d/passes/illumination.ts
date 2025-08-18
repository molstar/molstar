/**
 * Copyright (c) 2024-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CopyRenderable, QuadSchema, QuadValues, createCopyRenderable } from '../../mol-gl/compute/util';
import { DefineSpec, TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { Texture } from '../../mol-gl/webgl/texture';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { ValueCell } from '../../mol-util';
import { isDebugMode, isTimingMode } from '../../mol-util/debug';
import { Renderer, RendererProps } from '../../mol-gl/renderer';
import { Camera, ICamera } from '../camera';
import { Scene } from '../../mol-gl/scene';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { ShaderCode } from '../../mol-gl/shader-code';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { compose_frag } from '../../mol-gl/shader/illumination/compose.frag';
import { Vec2 } from '../../mol-math/linear-algebra/3d/vec2';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Vec3 } from '../../mol-math/linear-algebra/3d/vec3';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { Color } from '../../mol-util/color/color';
import { AntialiasingPass, PostprocessingPass, PostprocessingProps } from './postprocessing';
import { DrawPass } from './draw';
import { MarkingPass, MarkingProps } from './marking';
import { Helper } from '../helper/helper';
import { DofPass } from './dof';
import { TracingParams, TracingPass } from './tracing';
import { JitterVectors, MultiSampleProps } from './multi-sample';
import { compose_frag as multiSample_compose_frag } from '../../mol-gl/shader/compose.frag';
import { clamp, lerp } from '../../mol-math/interpolate';
import { SsaoProps } from './ssao';
import { OutlinePass } from './outline';
import { BloomPass } from './bloom';

type Props = {
    transparentBackground: boolean;
    dpoitIterations: number;
    illumination: IlluminationProps;
    renderer: RendererProps;
    postprocessing: PostprocessingProps;
    marking: MarkingProps;
    multiSample: MultiSampleProps;
}

type RenderContext = {
    renderer: Renderer;
    camera: Camera;
    scene: Scene;
    helper: Helper;
}

export const IlluminationParams = {
    enabled: PD.Boolean(false),
    maxIterations: PD.Numeric(5, { min: 0, max: 16, step: 1 }, { description: 'Maximum number of tracing iterations. Final iteration count is 2^x.' }),
    denoise: PD.Boolean(true),
    denoiseThreshold: PD.Interval([0.15, 1], { min: 0, max: 4, step: 0.01 }, { description: 'Threshold for denoising. Automatically adjusted within given interval based on current iteration.' }),
    ignoreOutline: PD.Boolean(true, { description: 'Ignore outline in illumination pass where it is generally not needed for visual clarity. Useful when illumination is often toggled on/off.' }),
    ...TracingParams,
};
export type IlluminationProps = PD.Values<typeof IlluminationParams>

export class IlluminationPass {
    private readonly tracing: TracingPass;

    private readonly transparentTarget: RenderTarget;
    private readonly outputTarget: RenderTarget;

    readonly packedDepth: boolean;

    private readonly copyRenderable: CopyRenderable;
    private readonly composeRenderable: ComposeRenderable;

    private multiSampleComposeTarget: RenderTarget;
    private multiSampleHoldTarget: RenderTarget;
    private multiSampleAccumulateTarget: RenderTarget;
    private multiSampleCompose: MultiSampleComposeRenderable;

    private _iteration = 0;
    get iteration() { return this._iteration; }

    private _colorTarget: RenderTarget;
    get colorTarget() { return this._colorTarget; }

    private _supported = false;
    get supported() {
        return this._supported;
    }

    getMaxIterations(props: Props) {
        return Math.pow(2, props.illumination.maxIterations);
    }

    static isSupported(webgl: WebGLContext) {
        const { drawBuffers, textureFloat, colorBufferFloat, depthTexture } = webgl.extensions;
        if (!textureFloat || !colorBufferFloat || !depthTexture || !drawBuffers) {
            if (isDebugMode) {
                const missing: string[] = [];
                if (!textureFloat) missing.push('textureFloat');
                if (!colorBufferFloat) missing.push('colorBufferFloat');
                if (!depthTexture) missing.push('depthTexture');
                if (!drawBuffers) missing.push('drawBuffers');
                console.log(`Missing "${missing.join('", "')}" extensions required for "illumination"`);
            }
            return false;
        } else {
            return true;
        }
    }

    constructor(private readonly webgl: WebGLContext, private readonly drawPass: DrawPass) {
        if (!IlluminationPass.isSupported(webgl)) return;

        const { colorTarget } = drawPass;
        const width = colorTarget.getWidth();
        const height = colorTarget.getHeight();

        this.tracing = new TracingPass(webgl, this.drawPass);

        this.transparentTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'nearest');
        this.outputTarget = webgl.createRenderTarget(width, height, false, 'uint8', 'linear');

        this.copyRenderable = createCopyRenderable(webgl, this.transparentTarget.texture);

        this.composeRenderable = getComposeRenderable(webgl, this.tracing.accumulateTarget.texture, this.tracing.normalTextureOpaque, this.tracing.colorTextureOpaque, this.drawPass.depthTextureOpaque, this.drawPass.depthTargetTransparent.texture, this.drawPass.postprocessing.outline.target.texture, this.transparentTarget.texture, this.drawPass.postprocessing.ssao.ssaoDepthTexture, this.drawPass.postprocessing.ssao.ssaoDepthTransparentTexture, false);

        this.multiSampleComposeTarget = webgl.createRenderTarget(width, height, false, 'float32');
        this.multiSampleHoldTarget = webgl.createRenderTarget(width, height, false);
        this.multiSampleAccumulateTarget = webgl.createRenderTarget(width, height, false, 'float32');
        this.multiSampleCompose = getMultiSampleComposeRenderable(webgl, this.outputTarget.texture);

        this._supported = true;
    }

    private renderInput(renderer: Renderer, camera: ICamera, scene: Scene, props: Props) {
        if (isTimingMode) this.webgl.timer.mark('IlluminationPass.renderInput');
        const { gl, state } = this.webgl;

        const markingEnabled = MarkingPass.isEnabled(props.marking);
        const hasTransparent = scene.opacityAverage < 1;
        const hasMarking = markingEnabled && scene.markerAverage > 0;

        this.transparentTarget.bind();
        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);

        const outlineEnabled = PostprocessingPass.isTransparentOutlineEnabled(props.postprocessing) && !props.illumination.ignoreOutline;
        const dofEnabled = DofPass.isEnabled(props.postprocessing);
        const ssaoEnabled = PostprocessingPass.isTransparentSsaoEnabled(scene, props.postprocessing);

        if (outlineEnabled || dofEnabled || ssaoEnabled) {
            this.drawPass.depthTargetTransparent.bind();
            renderer.clearDepth(true);
        }

        if (hasTransparent) {
            if (this.drawPass.transparency === 'wboit') {
                this.drawPass.wboit.bind();
                renderer.renderWboitTransparent(scene.primitives, camera, this.drawPass.depthTextureOpaque);

                if (scene.volumes.renderables.length > 0) {
                    renderer.renderWboitTransparent(scene.volumes, camera, this.drawPass.depthTextureOpaque);
                }

                this.transparentTarget.bind();
                this.drawPass.wboit.render();
            } else if (this.drawPass.transparency === 'dpoit') {
                const dpoitTextures = this.drawPass.dpoit.bind();
                renderer.renderDpoitTransparent(scene.primitives, camera, this.drawPass.depthTextureOpaque, dpoitTextures);

                for (let i = 0, iterations = props.dpoitIterations; i < iterations; i++) {
                    if (isTimingMode) this.webgl.timer.mark('DpoitPass.layer');
                    const dpoitTextures = this.drawPass.dpoit.bindDualDepthPeeling();
                    renderer.renderDpoitTransparent(scene.primitives, camera, this.drawPass.depthTextureOpaque, dpoitTextures);

                    if (iterations > 1) {
                        this.transparentTarget.bind();
                        this.drawPass.dpoit.renderBlendBack();
                    }
                    if (isTimingMode) this.webgl.timer.markEnd('DpoitPass.layer');
                }

                // evaluate dpoit
                this.transparentTarget.bind();
                this.drawPass.dpoit.render();

                if (scene.volumes.renderables.length > 0) {
                    renderer.renderVolume(scene.volumes, camera, this.drawPass.depthTextureOpaque);
                }
            } else {
                this.transparentTarget.bind();
                this.drawPass.depthTextureOpaque.attachFramebuffer(this.transparentTarget.framebuffer, 'depth');
                renderer.renderBlendedTransparent(scene.primitives, camera);
                this.drawPass.depthTextureOpaque.detachFramebuffer(this.transparentTarget.framebuffer, 'depth');

                if (scene.volumes.renderables.length > 0) {
                    renderer.renderVolume(scene.volumes, camera, this.drawPass.depthTextureOpaque);
                }
            }

            if (outlineEnabled || dofEnabled || ssaoEnabled) {
                this.drawPass.depthTargetTransparent.bind();
                if (scene.opacityAverage < 1) {
                    renderer.renderDepthTransparent(scene.primitives, camera, this.drawPass.depthTextureOpaque);
                }
            }

            if (ssaoEnabled) {
                this.drawPass.postprocessing.ssao.update(camera, scene, props.postprocessing.occlusion.params as SsaoProps, true);
                this.drawPass.postprocessing.ssao.render(camera);
            }
        }

        //

        if (hasMarking) {
            const markingDepthTest = props.marking.ghostEdgeStrength < 1;
            if (markingDepthTest && scene.markerAverage !== 1) {
                this.drawPass.marking.depthTarget.bind();
                renderer.clear(false, true);
                renderer.renderMarkingDepth(scene.primitives, camera);
            }

            this.drawPass.marking.maskTarget.bind();
            renderer.clear(false, true);
            renderer.renderMarkingMask(scene.primitives, camera, markingDepthTest ? this.drawPass.marking.depthTarget.texture : null);

            this.drawPass.marking.update(props.marking);
            this.drawPass.marking.render(camera.viewport, this.transparentTarget);
        }

        //

        this.tracing.composeTarget.bind();
        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);
        if (isTimingMode) this.webgl.timer.markEnd('IlluminationPass.renderInput');
    }

    shouldRender(props: Props) {
        return this._supported && props.illumination.enabled && this._iteration < this.getMaxIterations(props);
    }

    setSize(width: number, height: number) {
        if (!this._supported) return;

        const w = this.outputTarget.getWidth();
        const h = this.outputTarget.getHeight();

        if (width !== w || height !== h) {
            this.tracing.setSize(width, height);

            this.transparentTarget.setSize(width, height);
            this.outputTarget.setSize(width, height);

            ValueCell.update(this.copyRenderable.values.uTexSize, Vec2.set(this.copyRenderable.values.uTexSize.ref.value, width, height));
            ValueCell.update(this.composeRenderable.values.uTexSize, Vec2.set(this.composeRenderable.values.uTexSize.ref.value, width, height));

            this.multiSampleComposeTarget.setSize(width, height);
            this.multiSampleHoldTarget.setSize(width, height);
            this.multiSampleAccumulateTarget.setSize(width, height);
            ValueCell.update(this.multiSampleCompose.values.uTexSize, Vec2.set(this.multiSampleCompose.values.uTexSize.ref.value, width, height));
        }

        this.drawPass.setSize(width, height);
    }

    reset() {
        if (!this._supported) return;

        this.tracing.reset();
        this._iteration = 0;
        this.prevSampleIndex = -1;
    }

    restart(clearAdjustedProps = false) {
        if (!this._supported) return;

        this.tracing.restart(clearAdjustedProps);
        this._iteration = 0;
        this.prevSampleIndex = -1;
    }

    private renderInternal(ctx: RenderContext, props: Props, toDrawingBuffer: boolean, forceRenderInput: boolean) {
        if (!this.shouldRender(props)) return;

        if (isTimingMode) {
            this.webgl.timer.mark('IlluminationPass.render', {
                note: `iteration ${this._iteration + 1} of ${this.getMaxIterations(props)}`
            });
        }
        this.tracing.render(ctx, props.transparentBackground, props.illumination, this._iteration, forceRenderInput);

        const { renderer, camera, scene, helper } = ctx;
        const { gl, state } = this.webgl;
        const { x, y, width, height } = camera.viewport;

        if (this._iteration === 0 || forceRenderInput) {
            // render color & depth
            renderer.setTransparentBackground(props.transparentBackground);
            renderer.setDrawingBufferSize(this.tracing.composeTarget.getWidth(), this.tracing.composeTarget.getHeight());
            renderer.setPixelRatio(this.webgl.pixelRatio);
            renderer.setViewport(x, y, width, height);
            renderer.update(camera, scene);
            this.renderInput(renderer, camera, scene, props);
        }

        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.disable(gl.CULL_FACE);
        state.depthMask(false);
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

        const antialiasingEnabled = AntialiasingPass.isEnabled(props.postprocessing);
        const outlinesEnabled = OutlinePass.isEnabled(props.postprocessing) && !props.illumination.ignoreOutline;
        const occlusionEnabled = PostprocessingPass.isTransparentSsaoEnabled(scene, props.postprocessing);
        const bloomEnabled = BloomPass.isEnabled(props.postprocessing);
        const dofEnabled = DofPass.isEnabled(props.postprocessing);

        const markingEnabled = MarkingPass.isEnabled(props.marking);
        const hasTransparent = scene.opacityAverage < 1;
        const hasMarking = markingEnabled && scene.markerAverage > 0;

        let needsUpdateCompose = false;

        if (this.composeRenderable.values.dOutlineEnable.ref.value !== outlinesEnabled) {
            needsUpdateCompose = true;
            ValueCell.update(this.composeRenderable.values.dOutlineEnable, outlinesEnabled);
        }

        if (outlinesEnabled && props.postprocessing.outline.name === 'on') {
            const { transparentOutline, outlineScale } = this.drawPass.postprocessing.outline.update(camera, props.postprocessing.outline.params, this.drawPass.depthTargetTransparent.texture, this.drawPass.depthTextureOpaque);
            this.drawPass.postprocessing.outline.render();

            ValueCell.update(this.composeRenderable.values.uOutlineColor, Color.toVec3Normalized(this.composeRenderable.values.uOutlineColor.ref.value, props.postprocessing.outline.params.color));

            if (this.composeRenderable.values.dOutlineScale.ref.value !== outlineScale) {
                needsUpdateCompose = true;
                ValueCell.update(this.composeRenderable.values.dOutlineScale, outlineScale);
            }
            if (this.composeRenderable.values.dTransparentOutline.ref.value !== transparentOutline) {
                needsUpdateCompose = true;
                ValueCell.update(this.composeRenderable.values.dTransparentOutline, transparentOutline);
            }
        }

        if (this.composeRenderable.values.dOcclusionEnable.ref.value !== occlusionEnabled) {
            needsUpdateCompose = true;
            ValueCell.update(this.composeRenderable.values.dOcclusionEnable, occlusionEnabled);
        }

        if (occlusionEnabled && props.postprocessing.occlusion.name === 'on') {
            ValueCell.update(this.composeRenderable.values.uOcclusionColor, Color.toVec3Normalized(this.composeRenderable.values.uOcclusionColor.ref.value, props.postprocessing.occlusion.params.color));
        }

        const blendTransparency = hasTransparent || hasMarking;
        if (this.composeRenderable.values.dBlendTransparency.ref.value !== blendTransparency) {
            needsUpdateCompose = true;
            ValueCell.update(this.composeRenderable.values.dBlendTransparency, blendTransparency);
        }

        ValueCell.updateIfChanged(this.composeRenderable.values.uNear, camera.near);
        ValueCell.updateIfChanged(this.composeRenderable.values.uFar, camera.far);
        ValueCell.updateIfChanged(this.composeRenderable.values.uFogFar, camera.fogFar);
        ValueCell.updateIfChanged(this.composeRenderable.values.uFogNear, camera.fogNear);
        ValueCell.update(this.composeRenderable.values.uFogColor, Color.toVec3Normalized(this.composeRenderable.values.uFogColor.ref.value, renderer.props.backgroundColor));
        if (this.composeRenderable.values.dOrthographic.ref.value !== orthographic) {
            ValueCell.update(this.composeRenderable.values.dOrthographic, orthographic);
            needsUpdateCompose = true;
        }

        // background

        const _toDrawingBuffer = toDrawingBuffer && !antialiasingEnabled && !dofEnabled;
        if (_toDrawingBuffer) {
            this.webgl.bindDrawingBuffer();
        } else {
            this.tracing.composeTarget.bind();
        }
        this._colorTarget = this.tracing.composeTarget;

        this.drawPass.postprocessing.background.update(camera, props.postprocessing.background);
        this.drawPass.postprocessing.background.clear(props.postprocessing.background, props.transparentBackground, renderer.props.backgroundColor);
        this.drawPass.postprocessing.background.render(props.postprocessing.background);

        // compose

        ValueCell.updateIfChanged(this.composeRenderable.values.uTransparentBackground, props.transparentBackground || this.drawPass.postprocessing.background.isEnabled(props.postprocessing));
        if (this.composeRenderable.values.dDenoise.ref.value !== props.illumination.denoise) {
            ValueCell.update(this.composeRenderable.values.dDenoise, props.illumination.denoise);
            needsUpdateCompose = true;
        }
        const denoiseThreshold = props.multiSample.mode === 'on'
            ? props.illumination.denoiseThreshold[0]
            : lerp(props.illumination.denoiseThreshold[1], props.illumination.denoiseThreshold[0], clamp(this.iteration / (this.getMaxIterations(props) / 2), 0, 1));
        ValueCell.updateIfChanged(this.composeRenderable.values.uDenoiseThreshold, denoiseThreshold);
        if (needsUpdateCompose) this.composeRenderable.update();
        this.composeRenderable.render();

        //

        renderer.setDrawingBufferSize(this.tracing.composeTarget.getWidth(), this.tracing.composeTarget.getHeight());
        renderer.setPixelRatio(this.webgl.pixelRatio);
        renderer.setViewport(x, y, width, height);
        renderer.update(camera, scene);

        if (helper.debug.isEnabled) {
            helper.debug.syncVisibility();
            renderer.renderBlended(helper.debug.scene, camera);
        }
        if (helper.handle.isEnabled) {
            renderer.renderBlended(helper.handle.scene, camera);
        }
        if (helper.camera.isEnabled) {
            helper.camera.update(camera);
            renderer.update(helper.camera.camera, helper.camera.scene);
            renderer.renderBlended(helper.camera.scene, helper.camera.camera);
        }

        //

        let targetIsDrawingbuffer = false;
        let swapTarget = this.outputTarget;

        if (antialiasingEnabled) {
            const _toDrawingBuffer = toDrawingBuffer && !dofEnabled;
            this.drawPass.antialiasing.render(camera, this.tracing.composeTarget.texture, _toDrawingBuffer ? true : this.outputTarget, props.postprocessing);

            if (_toDrawingBuffer) {
                targetIsDrawingbuffer = true;
            } else {
                this._colorTarget = this.outputTarget;
                swapTarget = this.tracing.composeTarget;
            }
        }

        if (bloomEnabled && props.postprocessing.bloom.name === 'on') {
            const _toDrawingBuffer = (toDrawingBuffer && !dofEnabled) || targetIsDrawingbuffer;
            this.drawPass.bloom.update(this.tracing.colorTextureOpaque, this.tracing.normalTextureOpaque, this.drawPass.depthTextureOpaque, props.postprocessing.bloom.params);
            this.drawPass.bloom.render(camera.viewport, _toDrawingBuffer ? undefined : this._colorTarget);
        }

        if (dofEnabled && props.postprocessing.dof.name === 'on') {
            const _toDrawingBuffer = toDrawingBuffer || targetIsDrawingbuffer;
            this.drawPass.dof.update(camera, this._colorTarget.texture, this.drawPass.depthTextureOpaque, this.drawPass.depthTargetTransparent.texture, props.postprocessing.dof.params, scene.boundingSphereVisible);
            this.drawPass.dof.render(camera.viewport, _toDrawingBuffer ? undefined : swapTarget);

            if (!_toDrawingBuffer) {
                this._colorTarget = swapTarget;
            }
        }

        this._iteration += 1;
        if (isTimingMode) this.webgl.timer.markEnd('IlluminationPass.render');

        this.webgl.gl.flush();
    }

    private prevSampleIndex = -1;

    private renderMultiSample(ctx: RenderContext, props: Props, toDrawingBuffer: boolean) {
        const { camera } = ctx;
        const { multiSampleCompose, multiSampleComposeTarget, multiSampleHoldTarget, webgl } = this;
        const { gl, state } = webgl;

        // based on the Multisample Anti-Aliasing Render Pass
        // contributed to three.js by bhouston / http://clara.io/
        //
        // This manual approach to MSAA re-renders the scene once for
        // each sample with camera jitter and accumulates the results.
        const offsetList = JitterVectors[Math.max(0, Math.min(props.multiSample.sampleLevel, 5))];

        const maxIterations = this.getMaxIterations(props);
        const iteration = Math.min(this._iteration, maxIterations);

        const sampleIndex = Math.floor((iteration / maxIterations) * offsetList.length);

        if (isTimingMode) {
            webgl.timer.mark('IlluminationPass.renderMultiSample', {
                note: `sampleIndex ${sampleIndex + 1} of ${offsetList.length}`
            });
        }

        const { x, y, width, height } = camera.viewport;
        const sampleWeight = 1.0 / maxIterations;

        if (iteration === 0) {
            this.renderInternal(ctx, props, false, true);
            ValueCell.update(multiSampleCompose.values.uWeight, 1.0);
            ValueCell.update(multiSampleCompose.values.tColor, this._colorTarget.texture);
            multiSampleCompose.update();

            multiSampleHoldTarget.bind();
            state.disable(gl.BLEND);
            state.disable(gl.DEPTH_TEST);
            state.depthMask(false);
            state.viewport(x, y, width, height);
            state.scissor(x, y, width, height);
            multiSampleCompose.render();
        } else {
            camera.viewOffset.enabled = true;
            ValueCell.update(multiSampleCompose.values.tColor, this._colorTarget.texture);
            ValueCell.update(multiSampleCompose.values.uWeight, sampleWeight);
            multiSampleCompose.update();

            // render the scene multiple times, each slightly jitter offset
            // from the last and accumulate the results.
            const offset = offsetList[sampleIndex];
            Camera.setViewOffset(camera.viewOffset, width, height, offset[0], offset[1], width, height);
            camera.update();

            // render scene
            this.renderInternal(ctx, props, false, this.prevSampleIndex !== sampleIndex);

            // compose rendered scene with compose target
            multiSampleComposeTarget.bind();
            state.enable(gl.BLEND);
            state.blendEquationSeparate(gl.FUNC_ADD, gl.FUNC_ADD);
            state.blendFuncSeparate(gl.ONE, gl.ONE, gl.ONE, gl.ONE);
            state.disable(gl.DEPTH_TEST);
            state.depthMask(false);
            state.viewport(x, y, width, height);
            state.scissor(x, y, width, height);
            if (iteration === 1) {
                state.clearColor(0, 0, 0, 0);
                gl.clear(gl.COLOR_BUFFER_BIT);
            }
            multiSampleCompose.render();
        }

        this.prevSampleIndex = sampleIndex;

        if (toDrawingBuffer) {
            this.webgl.bindDrawingBuffer();
        } else {
            this.multiSampleAccumulateTarget.bind();
        }
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        const accumulationWeight = iteration * sampleWeight;
        if (accumulationWeight > 0) {
            ValueCell.update(multiSampleCompose.values.uWeight, 1.0);
            ValueCell.update(multiSampleCompose.values.tColor, multiSampleComposeTarget.texture);
            multiSampleCompose.update();
            state.disable(gl.BLEND);
            multiSampleCompose.render();
        }
        if (accumulationWeight < 1.0) {
            ValueCell.update(multiSampleCompose.values.uWeight, 1.0 - accumulationWeight);
            ValueCell.update(multiSampleCompose.values.tColor, multiSampleHoldTarget.texture);
            multiSampleCompose.update();
            if (accumulationWeight === 0) state.disable(gl.BLEND);
            else state.enable(gl.BLEND);
            multiSampleCompose.render();
        }

        if (!toDrawingBuffer) {
            state.disable(gl.BLEND);
            this.colorTarget.bind();
            if (this.copyRenderable.values.tColor.ref.value !== this.multiSampleAccumulateTarget.texture) {
                ValueCell.update(this.copyRenderable.values.tColor, this.multiSampleAccumulateTarget.texture);
                this.copyRenderable.update();
            }
            this.copyRenderable.render();
        }

        camera.viewOffset.enabled = false;
        camera.update();
        if (isTimingMode) webgl.timer.markEnd('IlluminationPass.renderMultiSample');
    }

    render(ctx: RenderContext, props: Props, toDrawingBuffer: boolean) {
        if (!this._supported) return;

        if (props.multiSample.mode === 'on') {
            this.renderMultiSample(ctx, props, toDrawingBuffer);
        } else {
            this.renderInternal(ctx, props, toDrawingBuffer, false);
        }
    }
}

//

const ComposeSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tNormal: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tShaded: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tTransparentColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    dBlendTransparency: DefineSpec('boolean'),
    tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tSsaoDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tOutlines: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),

    dDenoise: DefineSpec('boolean'),
    uDenoiseThreshold: UniformSpec('f'),

    dOrthographic: DefineSpec('number'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
    uFogColor: UniformSpec('v3'),
    uOutlineColor: UniformSpec('v3'),
    uOcclusionColor: UniformSpec('v3'),
    uTransparentBackground: UniformSpec('b'),

    dOcclusionEnable: DefineSpec('boolean'),
    dOutlineEnable: DefineSpec('boolean'),
    dOutlineScale: DefineSpec('number'),
    dTransparentOutline: DefineSpec('boolean'),
};
const ComposeShaderCode = ShaderCode('compose', quad_vert, compose_frag);
type ComposeRenderable = ComputeRenderable<Values<typeof ComposeSchema>>

function getComposeRenderable(ctx: WebGLContext, colorTexture: Texture, normalTexture: Texture, shadedTexture: Texture, depthTextureOpaque: Texture, depthTextureTransparent: Texture, outlinesTexture: Texture, transparentColorTexture: Texture, ssaoDepthOpaqueTexture: Texture, ssaoDepthTransparentTexture: Texture, transparentOutline: boolean): ComposeRenderable {
    const values: Values<typeof ComposeSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        tNormal: ValueCell.create(normalTexture),
        tShaded: ValueCell.create(shadedTexture),
        tTransparentColor: ValueCell.create(transparentColorTexture),
        dBlendTransparency: ValueCell.create(true),
        tSsaoDepth: ValueCell.create(ssaoDepthOpaqueTexture),
        tSsaoDepthTransparent: ValueCell.create(ssaoDepthTransparentTexture),
        tDepthOpaque: ValueCell.create(depthTextureOpaque),
        tDepthTransparent: ValueCell.create(depthTextureTransparent),
        tOutlines: ValueCell.create(outlinesTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),

        dDenoise: ValueCell.create(true),
        uDenoiseThreshold: ValueCell.create(0.1),

        dOrthographic: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uFogNear: ValueCell.create(10000),
        uFogFar: ValueCell.create(10000),
        uFogColor: ValueCell.create(Vec3.create(1, 1, 1)),
        uOutlineColor: ValueCell.create(Vec3.create(0, 0, 0)),
        uOcclusionColor: ValueCell.create(Vec3.create(0, 0, 0)),
        uTransparentBackground: ValueCell.create(false),

        dOcclusionEnable: ValueCell.create(false),
        dOutlineEnable: ValueCell.create(false),
        dOutlineScale: ValueCell.create(1),
        dTransparentOutline: ValueCell.create(transparentOutline),
    };

    const schema = { ...ComposeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', ComposeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const MultiSampleComposeSchema = {
    ...QuadSchema,
    tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    uTexSize: UniformSpec('v2'),
    uWeight: UniformSpec('f'),
};
const MultiSampleComposeShaderCode = ShaderCode('compose', quad_vert, multiSample_compose_frag);
type MultiSampleComposeRenderable = ComputeRenderable<Values<typeof MultiSampleComposeSchema>>

function getMultiSampleComposeRenderable(ctx: WebGLContext, colorTexture: Texture): MultiSampleComposeRenderable {
    const values: Values<typeof MultiSampleComposeSchema> = {
        ...QuadValues,
        tColor: ValueCell.create(colorTexture),
        uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),
        uWeight: ValueCell.create(1.0),
    };

    const schema = { ...MultiSampleComposeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', MultiSampleComposeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}
