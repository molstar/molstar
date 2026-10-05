/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { type CopyRenderable, createCopyRenderable, getSharedCopyRenderable, QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { WebGLContext } from '@molstar/graphics/gl/webgl/context';
import { createComputeRenderItem } from '@molstar/graphics/gl/webgl/render-item';
import type { Texture } from '@molstar/graphics/gl/webgl/texture';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { quad_vert } from '@molstar/graphics/gl/shader/quad.vert';
import { overlay_frag } from '@molstar/graphics/gl/shader/marking/overlay.frag';
import type { Viewport } from '../camera/util.js';
import type { RenderTarget } from '@molstar/graphics/gl/webgl/render-target';
import { Color } from '@molstar/core/util/color';
import { edge_frag } from '@molstar/graphics/gl/shader/marking/edge.frag';
import { dim_frag } from '@molstar/graphics/gl/shader/marking/dim.frag';
import { isTimingMode } from '@molstar/core/util/debug';
import { AntialiasingPass, type PostprocessingProps } from './postprocessing.js';
import { Camera, type ICamera } from '../camera.js';
import { StereoCamera } from '../camera/stereo.js';
import type { RendererProps } from '@molstar/graphics/gl/renderer';
import { Scene } from '@molstar/graphics/gl/scene';
import { clearJitter, setJitter } from './jitter.js';
import type { RenderContext as BaseRenderContext } from '../util.js';

export const MarkingParams = {
    enabled: PD.Boolean(true),
    highlightEdgeColor: PD.Color(Color.darken(Color.fromNormalizedRgb(1.0, 0.4, 0.6), 1.0)),
    selectEdgeColor: PD.Color(Color.darken(Color.fromNormalizedRgb(0.2, 1.0, 0.1), 1.0)),
    edgeScale: PD.Numeric(1, { min: 1, max: 3, step: 0.1 }, { description: 'Thickness of the edge.' }),
    highlightEdgeStrength: PD.Numeric(1.0, { min: 0, max: 1, step: 0.1 }),
    selectEdgeStrength: PD.Numeric(1.0, { min: 0, max: 1, step: 0.1 }),
    ghostEdgeStrength: PD.Numeric(0.3, { min: 0, max: 1, step: 0.1 }, { description: 'Opacity of the hidden edges that are covered by other geometry. When set to 1 and not dimming, one less geometry render pass is done, but hidden parts are then also tinted.' }),
    innerEdgeFactor: PD.Numeric(1.5, { min: 0, max: 3, step: 0.1 }, { description: 'Factor to multiply the inner edge color with - for added contrast.' }),
};
export type MarkingProps = PD.Values<typeof MarkingParams>

export const SingleSample: number[][] = [[0, 0]];

/** Whether marked objects are tinted and unmarked ones dimmed in the material color instead of by the marking pass. */
export function isMaterialColorMarker(rendererProps: Pick<RendererProps, 'colorMarker'>, markingProps: MarkingProps) {
    return rendererProps.colorMarker && !markingProps.enabled;
}

type Props = {
    marking: MarkingProps;
    postprocessing: PostprocessingProps;
    transparentBackground: boolean;
}

export type MarkingPresentOptions = {
    /** image to blend the marking over */
    base: RenderTarget
    toDrawingBuffer: boolean
    /** jitter offsets the layer is accumulated from */
    offsets: number[][]
    /** start the layer over instead of adding to it */
    restart: boolean
    /** maximum number of samples to add */
    samples: number
    /** how the base image was shaded, kept on tinted and dimmed objects */
    shading: MarkingShading | null
}

export type MarkingShading =
    | { name: 'ssao', ssao: Texture }
    /** the (denoised) base image relative to the direct `shaded` color gives the path traced occlusion and shadows */
    | { name: 'traced', shaded: Texture }

type RenderContext = BaseRenderContext<ICamera | StereoCamera>

export class MarkingPass {
    static isEnabled(props: MarkingProps) {
        return props.enabled;
    }

    static hasMarking(scene: Scene, props: { marking: MarkingProps }) {
        return MarkingPass.isEnabled(props.marking) && scene.markerAverage > 0;
    }

    private readonly depthTarget: RenderTarget;
    private readonly maskTarget: RenderTarget;
    private readonly edgesTarget: RenderTarget;
    /** premultiplied marking, the running average of the samples taken so far */
    private readonly layerTarget: RenderTarget;

    private readonly edge: EdgeRenderable;
    private readonly overlay: OverlayRenderable;
    private readonly composite: CopyRenderable;

    // the mask is antialiased with the same method as the rest of the image, created on first use
    private antialiasing: AntialiasingPass | null = null;

    /** coverage of the visible unmarked geometry, created on first use */
    private dim: { target: RenderTarget, aaTarget: RenderTarget, renderable: DimRenderable } | null = null;

    /** jitter offsets the layer is accumulated from */
    private offsets = SingleSample;
    /** number of offsets the layer holds samples for, 0 if it is not valid */
    private sampleCount = 0;

    /** image of the last render to the drawing buffer, before marking was blended over it */
    private base: RenderTarget | undefined = undefined;
    /** shading of `base` */
    private baseShading: MarkingShading | null = null;

    constructor(private webgl: WebGLContext, width: number, height: number) {
        const { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } = webgl.extensions;
        this.depthTarget = webgl.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'nearest');
        // linear so that it can be used as antialiasing input
        this.maskTarget = webgl.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'linear');
        this.edgesTarget = webgl.createRenderTarget(width, height);
        const layerType = colorBufferHalfFloat && textureHalfFloat ? 'fp16' :
            colorBufferFloat && textureFloat ? 'float32' : 'uint8';
        this.layerTarget = webgl.createRenderTarget(width, height, 'none', layerType);

        this.edge = getEdgeRenderable(webgl, this.maskTarget.texture);
        this.overlay = getOverlayRenderable(webgl, this.edgesTarget.texture, this.maskTarget.texture);
        this.composite = createCopyRenderable(webgl, this.layerTarget.texture);
    }

    getByteCount() {
        const dimByteCount = this.dim ? this.dim.target.getByteCount() + this.dim.aaTarget.getByteCount() : 0;
        return this.depthTarget.getByteCount() + this.maskTarget.getByteCount() + this.edgesTarget.getByteCount() + this.layerTarget.getByteCount() + (this.antialiasing?.getByteCount() ?? 0) + dimByteCount;
    }

    private setEdgeState(viewport: Viewport) {
        const { gl, state } = this.webgl;

        state.enable(gl.SCISSOR_TEST);
        state.enable(gl.BLEND);
        state.blendFunc(gl.ONE, gl.ONE);
        state.blendEquation(gl.FUNC_ADD);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);

        const { x, y, width, height } = viewport;
        state.viewport(x, y, width, height);
        state.scissor(x, y, width, height);

        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);
    }

    private setViewport(viewport: Viewport) {
        const { x, y, width, height } = viewport;
        this.webgl.state.viewport(x, y, width, height);
        this.webgl.state.scissor(x, y, width, height);
    }

    setSize(width: number, height: number) {
        const w = this.depthTarget.getWidth();
        const h = this.depthTarget.getHeight();

        if (width !== w || height !== h) {
            this.depthTarget.setSize(width, height);
            this.maskTarget.setSize(width, height);
            this.edgesTarget.setSize(width, height);
            this.layerTarget.setSize(width, height);
            this.antialiasing?.setSize(width, height);
            if (this.dim) {
                this.dim.target.setSize(width, height);
                this.dim.aaTarget.setSize(width, height);
                ValueCell.update(this.dim.renderable.values.uTexSizeInv, Vec2.set(this.dim.renderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            }
            this.sampleCount = 0;
            this.invalidate();

            ValueCell.update(this.edge.values.uTexSizeInv, Vec2.set(this.edge.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.overlay.values.uTexSizeInv, Vec2.set(this.overlay.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.composite.values.uTexSize, Vec2.set(this.composite.values.uTexSize.ref.value, width, height));
        }
    }

    private update(props: Props, rendererProps: RendererProps, dim: boolean, depthTest: boolean, shading: MarkingShading | null, base: Texture) {
        const { highlightEdgeColor, selectEdgeColor, edgeScale, innerEdgeFactor, ghostEdgeStrength, highlightEdgeStrength, selectEdgeStrength } = props.marking;

        const { values: edgeValues } = this.edge;
        ValueCell.updateIfChanged(edgeValues.uEdgeScale, Math.max(1, edgeScale * this.webgl.pixelRatio));

        const { values: overlayValues } = this.overlay;
        ValueCell.update(overlayValues.uHighlightEdgeColor, Color.toVec3Normalized(overlayValues.uHighlightEdgeColor.ref.value, highlightEdgeColor));
        ValueCell.update(overlayValues.uSelectEdgeColor, Color.toVec3Normalized(overlayValues.uSelectEdgeColor.ref.value, selectEdgeColor));
        ValueCell.updateIfChanged(overlayValues.uInnerEdgeFactor, innerEdgeFactor);
        ValueCell.updateIfChanged(overlayValues.uGhostEdgeStrength, ghostEdgeStrength);
        ValueCell.updateIfChanged(overlayValues.uDepthTest, depthTest);
        ValueCell.updateIfChanged(overlayValues.uHighlightEdgeStrength, highlightEdgeStrength);
        ValueCell.updateIfChanged(overlayValues.uSelectEdgeStrength, selectEdgeStrength);

        const { colorMarker, highlightColor, selectColor, highlightStrength, selectStrength, dimColor, dimStrength, backgroundColor } = rendererProps;
        const fill = colorMarker && !isMaterialColorMarker(rendererProps, props.marking);
        ValueCell.update(overlayValues.uHighlightFillColor, Color.toVec3Normalized(overlayValues.uHighlightFillColor.ref.value, highlightColor));
        ValueCell.update(overlayValues.uSelectFillColor, Color.toVec3Normalized(overlayValues.uSelectFillColor.ref.value, selectColor));
        ValueCell.updateIfChanged(overlayValues.uHighlightFillStrength, fill ? highlightStrength : 0);
        ValueCell.updateIfChanged(overlayValues.uSelectFillStrength, fill ? selectStrength : 0);
        ValueCell.update(overlayValues.uDimColor, Color.toVec3Normalized(overlayValues.uDimColor.ref.value, dimColor));
        ValueCell.updateIfChanged(overlayValues.uDimStrength, dim ? dimStrength : 0);

        const occlusionProps = props.postprocessing.occlusion;
        const shaded = dim || (fill && (highlightStrength > 0 || selectStrength > 0));
        const markingShading = !shaded || !shading ? 'off'
            : shading.name === 'traced' ? 'traced'
                : occlusionProps.name === 'on' ? 'ssao' : 'off';
        let needsUpdate = false;
        if (overlayValues.dMarkingShading.ref.value !== markingShading) {
            ValueCell.update(overlayValues.dMarkingShading, markingShading);
            needsUpdate = true;
        }
        if (shading?.name === 'ssao' && overlayValues.tSsaoDepth.ref.value !== shading.ssao) {
            ValueCell.update(overlayValues.tSsaoDepth, shading.ssao);
            needsUpdate = true;
        }
        if (shading?.name === 'traced' && (overlayValues.tBase.ref.value !== base || overlayValues.tShaded.ref.value !== shading.shaded)) {
            ValueCell.update(overlayValues.tBase, base);
            ValueCell.update(overlayValues.tShaded, shading.shaded);
            needsUpdate = true;
        }
        if (needsUpdate) this.overlay.update();

        if (markingShading === 'ssao' && occlusionProps.name === 'on') {
            ValueCell.update(overlayValues.uOcclusionColor, Color.toVec3Normalized(overlayValues.uOcclusionColor.ref.value, occlusionProps.params.color));
        }
        if (markingShading !== 'off') {
            ValueCell.update(overlayValues.uFogColor, Color.toVec3Normalized(overlayValues.uFogColor.ref.value, backgroundColor));
            ValueCell.updateIfChanged(overlayValues.uTransparentBackground, props.transparentBackground);
        }
    }

    /** Forgets the base image, needed when the last render did not go through `present`. */
    invalidate() {
        this.base = undefined;
        this.baseShading = null;
    }

    /** Whether the layer is still missing samples that can be added by `redraw`. */
    get needsMoreSamples() {
        return !!this.base && this.sampleCount > 0 && this.sampleCount < this.offsets.length;
    }

    /**
     * Adds up to `samples` jittered samples to the marking layer (starting over if `restart`)
     * and blends it over `base`, either into the drawing buffer (remembering `base` for redraws)
     * or in place.
     */
    present(ctx: RenderContext, props: Props, options: MarkingPresentOptions) {
        const { base, toDrawingBuffer, shading } = options;
        const { viewport } = ctx.camera;
        const hasMarking = this.updateLayer(ctx, props, options);

        if (toDrawingBuffer) this.copyToDrawingBuffer(base, viewport);
        else if (hasMarking) base.bind();
        if (hasMarking) this.compositeLayer(viewport);

        if (toDrawingBuffer && !(ctx.camera instanceof StereoCamera)) {
            this.base = base;
            this.baseShading = shading;
        } else {
            this.invalidate();
        }
    }

    /**
     * Redraws marking over the image of the last render to the drawing buffer, without
     * re-rendering the scene. Returns false if that image is not available.
     */
    redraw(ctx: RenderContext, props: Props, offsets: number[][], restart: boolean, samples: number): boolean {
        if (!this.base) return false;

        if (isTimingMode) this.webgl.timer.mark('MarkingPass.redraw');
        this.present(ctx, props, { base: this.base, toDrawingBuffer: true, offsets, restart, samples, shading: this.baseShading });
        this.webgl.gl.flush();
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.redraw');
        return true;
    }

    /** Returns false if there is no marking. */
    private updateLayer(ctx: RenderContext, props: Props, options: MarkingPresentOptions): boolean {
        const { base, offsets, restart, samples, shading } = options;
        const { renderer, camera, scene, frame } = ctx;
        if (camera instanceof StereoCamera || camera.disabled || !MarkingPass.hasMarking(scene, props)) {
            this.sampleCount = 0;
            return false;
        }
        if (restart || offsets !== this.offsets) {
            this.offsets = offsets;
            this.sampleCount = 0;
        }

        const { colorMarker, dimStrength } = renderer.props;
        // nothing to dim when everything is marked
        const dim = colorMarker && dimStrength > 0 && scene.markerAverage < 1;
        // dimming needs to know where marked objects are visible
        const depthTest = props.marking.ghostEdgeStrength < 1 || dim;

        const end = Math.min(this.sampleCount + samples, offsets.length);
        if (this.sampleCount >= end) {
            // traced shading comes from the base image, which keeps converging; the textures of a single sample are still there
            if (shading?.name === 'traced' && offsets.length === 1) {
                this.update(props, renderer.props, dim, depthTest, shading, base.texture);
                this.renderOverlay(ctx.camera.viewport, 1);
            }
            return true;
        }

        if (isTimingMode) this.webgl.timer.mark('MarkingPass.updateLayer');
        const jitter = offsets.length > 1 && camera instanceof Camera;
        const { x, y, width, height } = camera.viewport;

        renderer.setDrawingBufferSize(this.maskTarget.getWidth(), this.maskTarget.getHeight());
        renderer.setPixelRatio(this.webgl.pixelRatio);
        renderer.setViewport(x, y, width, height);
        this.update(props, renderer.props, dim, depthTest, shading, base.texture);

        for (; this.sampleCount < end; ++this.sampleCount) {
            const offset = offsets[this.sampleCount];
            if (jitter) setJitter(camera, offset);
            renderer.update(camera, scene, frame);

            if (depthTest) {
                this.depthTarget.bind();
                renderer.clear(false, true);
                // when everything is marked there is nothing to occlude with, but the cleared depth is still needed
                if (scene.markerAverage !== 1) renderer.renderMarkingDepth(scene.primitives, camera);
            }

            this.maskTarget.bind();
            renderer.clear(false, true);
            renderer.renderMarkingMask(scene.primitives, camera, depthTest ? this.depthTarget.texture : null);

            // the scene SSAO is computed without jitter, see `MultiSamplePass`
            ValueCell.update(this.overlay.values.uOcclusionOffset, Vec2.set(this.overlay.values.uOcclusionOffset.ref.value, jitter ? offset[0] / width : 0, jitter ? offset[1] / height : 0));
            this.renderSample(camera, props.postprocessing, 1 / (this.sampleCount + 1), dim);
        }

        if (jitter) clearJitter(camera);
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.updateLayer');
        return true;
    }

    private getAntialiasing(props: PostprocessingProps): AntialiasingPass | null {
        // no sharpening, it would break decoding the coverage from the mask
        if (!props.enabled || props.antialiasing.name === 'off') return null;

        if (!this.antialiasing) {
            // linear, as the edge pass samples the mask between texels
            this.antialiasing = new AntialiasingPass(this.webgl, this.maskTarget.getWidth(), this.maskTarget.getHeight(), 'linear');
        }
        return this.antialiasing;
    }

    /** Antialiases `input` into `output`, returns the texture to use. */
    private antialias(camera: ICamera, props: PostprocessingProps, input: RenderTarget, output: RenderTarget): Texture {
        const aa = this.getAntialiasing(props);
        return aa?.renderAntialiasingOnly(camera, input.texture, output, props) ? output.texture : input.texture;
    }

    private getDim() {
        if (!this.dim) {
            const width = this.maskTarget.getWidth();
            const height = this.maskTarget.getHeight();
            this.dim = {
                // linear so that it can be used as antialiasing input
                target: this.webgl.createRenderTarget(width, height, 'none', 'uint8', 'linear'),
                aaTarget: this.webgl.createRenderTarget(width, height, 'none'),
                renderable: getDimRenderable(this.webgl, this.depthTarget.texture, this.maskTarget.texture),
            };
        }
        return this.dim;
    }

    /** Renders the antialiased coverage of the visible unmarked geometry, from the depth and raw mask of the current sample. */
    private renderDim(camera: ICamera, postprocessingProps: PostprocessingProps): Texture {
        const { gl, state } = this.webgl;
        const { target, aaTarget, renderable } = this.getDim();
        const { values } = renderable;
        ValueCell.updateIfChanged(values.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
        ValueCell.updateIfChanged(values.uNear, camera.near);
        ValueCell.updateIfChanged(values.uFar, camera.far);
        ValueCell.updateIfChanged(values.uFogNear, camera.fogNear);
        ValueCell.updateIfChanged(values.uFogFar, camera.fogFar);

        target.bind();
        this.setViewport(camera.viewport);
        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
        state.clearColor(0, 0, 0, 0);
        gl.clear(gl.COLOR_BUFFER_BIT);
        renderable.render();

        return this.antialias(camera, postprocessingProps, target, aaTarget);
    }

    /** Blends the marking of the current mask into the layer, with `weight` for the new sample. */
    private renderSample(camera: ICamera, postprocessingProps: PostprocessingProps, weight: number, dim: boolean) {
        if (isTimingMode) this.webgl.timer.mark('MarkingPass.renderSample');

        const dimTexture = dim ? this.renderDim(camera, postprocessingProps) : this.maskTarget.texture;
        const aa = this.getAntialiasing(postprocessingProps);
        const maskTexture = aa ? this.antialias(camera, postprocessingProps, this.maskTarget, aa.target) : this.maskTarget.texture;

        if (this.edge.values.tMaskTexture.ref.value !== maskTexture) {
            ValueCell.update(this.edge.values.tMaskTexture, maskTexture);
            this.edge.update();
            ValueCell.update(this.overlay.values.tMaskTexture, maskTexture);
            this.overlay.update();
        }
        if (this.overlay.values.tDimTexture.ref.value !== dimTexture) {
            ValueCell.update(this.overlay.values.tDimTexture, dimTexture);
            this.overlay.update();
        }

        const viewport = camera.viewport;
        this.edgesTarget.bind();
        this.setEdgeState(viewport);
        this.edge.render();

        this.renderOverlay(viewport, weight);
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.renderSample');
    }

    /** Blends the overlay of the current sample's textures into the layer, with `weight` for the new sample. */
    private renderOverlay(viewport: Viewport, weight: number) {
        const { gl, state } = this.webgl;
        this.layerTarget.bind();
        this.setViewport(viewport);
        state.enable(gl.SCISSOR_TEST);
        state.enable(gl.BLEND);
        // running average: the first sample (weight 1) replaces the previous content
        state.blendColor(0, 0, 0, weight);
        state.blendFunc(gl.CONSTANT_ALPHA, gl.ONE_MINUS_CONSTANT_ALPHA);
        state.blendEquation(gl.FUNC_ADD);
        this.overlay.render();
    }

    /** Blends the layer over the currently bound framebuffer. */
    private compositeLayer(viewport: Viewport) {
        const { gl, state } = this.webgl;
        state.enable(gl.SCISSOR_TEST);
        state.enable(gl.BLEND);
        state.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
        state.blendEquation(gl.FUNC_ADD);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
        this.setViewport(viewport);
        this.composite.render();
    }

    private copyToDrawingBuffer(src: RenderTarget, viewport: Viewport) {
        const { gl, state } = this.webgl;
        const copy = getSharedCopyRenderable(this.webgl, src.texture);
        this.webgl.bindDrawingBuffer();
        state.enable(gl.SCISSOR_TEST);
        state.disable(gl.BLEND);
        state.disable(gl.DEPTH_TEST);
        state.depthMask(false);
        this.setViewport(viewport);
        copy.render();
    }
}

//

const EdgeSchema = {
    ...QuadSchema,
    tMaskTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uEdgeScale: UniformSpec('f'),
};
const EdgeShaderCode = ShaderCode('edge', quad_vert, edge_frag);
type EdgeRenderable = ComputeRenderable<Values<typeof EdgeSchema>>

function getEdgeRenderable(ctx: WebGLContext, maskTexture: Texture): EdgeRenderable {
    const width = maskTexture.getWidth();
    const height = maskTexture.getHeight();

    const values: Values<typeof EdgeSchema> = {
        ...QuadValues,
        tMaskTexture: ValueCell.create(maskTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uEdgeScale: ValueCell.create(1),
    };

    const schema = { ...EdgeSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', EdgeShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const DimSchema = {
    ...QuadSchema,
    tDepthTexture: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tMaskTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uIsOrtho: UniformSpec('f'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
};
const DimShaderCode = ShaderCode('dim', quad_vert, dim_frag);
type DimRenderable = ComputeRenderable<Values<typeof DimSchema>>

function getDimRenderable(ctx: WebGLContext, depthTexture: Texture, maskTexture: Texture): DimRenderable {
    const width = maskTexture.getWidth();
    const height = maskTexture.getHeight();

    const values: Values<typeof DimSchema> = {
        ...QuadValues,
        tDepthTexture: ValueCell.create(depthTexture),
        tMaskTexture: ValueCell.create(maskTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uIsOrtho: ValueCell.create(0),
        uNear: ValueCell.create(1),
        uFar: ValueCell.create(10000),
        uFogNear: ValueCell.create(1),
        uFogFar: ValueCell.create(10000),
    };

    const schema = { ...DimSchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', DimShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}

//

const OverlaySchema = {
    ...QuadSchema,
    tEdgeTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tMaskTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tDimTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uHighlightEdgeColor: UniformSpec('v3'),
    uSelectEdgeColor: UniformSpec('v3'),
    uHighlightEdgeStrength: UniformSpec('f'),
    uSelectEdgeStrength: UniformSpec('f'),
    uGhostEdgeStrength: UniformSpec('f'),
    uDepthTest: UniformSpec('b'),
    uInnerEdgeFactor: UniformSpec('f'),
    uHighlightFillColor: UniformSpec('v3'),
    uSelectFillColor: UniformSpec('v3'),
    uHighlightFillStrength: UniformSpec('f'),
    uSelectFillStrength: UniformSpec('f'),
    uDimColor: UniformSpec('v3'),
    uDimStrength: UniformSpec('f'),
    tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uOcclusionColor: UniformSpec('v3'),
    uFogColor: UniformSpec('v3'),
    uTransparentBackground: UniformSpec('b'),
    uOcclusionOffset: UniformSpec('v2'),
    tBase: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    tShaded: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
    dMarkingShading: DefineSpec('string', ['off', 'ssao', 'traced']),
};
const OverlayShaderCode = ShaderCode('overlay', quad_vert, overlay_frag);
type OverlayRenderable = ComputeRenderable<Values<typeof OverlaySchema>>

function getOverlayRenderable(ctx: WebGLContext, edgeTexture: Texture, maskTexture: Texture): OverlayRenderable {
    const width = edgeTexture.getWidth();
    const height = edgeTexture.getHeight();

    const values: Values<typeof OverlaySchema> = {
        ...QuadValues,
        tEdgeTexture: ValueCell.create(edgeTexture),
        tMaskTexture: ValueCell.create(maskTexture),
        // placeholder until dimming is used
        tDimTexture: ValueCell.create(maskTexture),
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uHighlightEdgeColor: ValueCell.create(Vec3()),
        uSelectEdgeColor: ValueCell.create(Vec3()),
        uHighlightEdgeStrength: ValueCell.create(1),
        uSelectEdgeStrength: ValueCell.create(1),
        uGhostEdgeStrength: ValueCell.create(0),
        uDepthTest: ValueCell.create(false),
        uInnerEdgeFactor: ValueCell.create(0),
        uHighlightFillColor: ValueCell.create(Vec3()),
        uSelectFillColor: ValueCell.create(Vec3()),
        uHighlightFillStrength: ValueCell.create(0),
        uSelectFillStrength: ValueCell.create(0),
        uDimColor: ValueCell.create(Vec3()),
        uDimStrength: ValueCell.create(0),
        // placeholders until marking shading is used
        tSsaoDepth: ValueCell.create(maskTexture),
        uOcclusionColor: ValueCell.create(Vec3()),
        uFogColor: ValueCell.create(Vec3()),
        uTransparentBackground: ValueCell.create(false),
        uOcclusionOffset: ValueCell.create(Vec2()),
        tBase: ValueCell.create(maskTexture),
        tShaded: ValueCell.create(maskTexture),
        dMarkingShading: ValueCell.create('off'),
    };

    const schema = { ...OverlaySchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', OverlayShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}