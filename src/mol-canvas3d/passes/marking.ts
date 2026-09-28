/**
 * Copyright (c) 2021-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { CopyRenderable, createCopyRenderable, getSharedCopyRenderable, QuadSchema, QuadValues } from '../../mol-gl/compute/util';
import { ComputeRenderable, createComputeRenderable } from '../../mol-gl/renderable';
import { TextureSpec, UniformSpec, Values } from '../../mol-gl/renderable/schema';
import { ShaderCode } from '../../mol-gl/shader-code';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { createComputeRenderItem } from '../../mol-gl/webgl/render-item';
import { Texture } from '../../mol-gl/webgl/texture';
import { Vec2, Vec3 } from '../../mol-math/linear-algebra';
import { ValueCell } from '../../mol-util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { quad_vert } from '../../mol-gl/shader/quad.vert';
import { overlay_frag } from '../../mol-gl/shader/marking/overlay.frag';
import { Viewport } from '../camera/util';
import { RenderTarget } from '../../mol-gl/webgl/render-target';
import { Color } from '../../mol-util/color';
import { edge_frag } from '../../mol-gl/shader/marking/edge.frag';
import { isTimingMode } from '../../mol-util/debug';
import { AntialiasingPass, PostprocessingProps } from './postprocessing';
import { Camera, ICamera } from '../camera';
import { StereoCamera } from '../camera/stereo';
import { RendererProps } from '../../mol-gl/renderer';
import { Scene } from '../../mol-gl/scene';
import { clearJitter, setJitter } from './jitter';
import { RenderContext as BaseRenderContext } from '../util';

export const MarkingParams = {
    enabled: PD.Boolean(true),
    highlightEdgeColor: PD.Color(Color.darken(Color.fromNormalizedRgb(1.0, 0.4, 0.6), 1.0)),
    selectEdgeColor: PD.Color(Color.darken(Color.fromNormalizedRgb(0.2, 1.0, 0.1), 1.0)),
    edgeScale: PD.Numeric(1, { min: 1, max: 3, step: 0.1 }, { description: 'Thickness of the edge.' }),
    highlightEdgeStrength: PD.Numeric(1.0, { min: 0, max: 1, step: 0.1 }),
    selectEdgeStrength: PD.Numeric(1.0, { min: 0, max: 1, step: 0.1 }),
    ghostEdgeStrength: PD.Numeric(0.3, { min: 0, max: 1, step: 0.1 }, { description: 'Opacity of the hidden edges that are covered by other geometry. When set to 1, one less geometry render pass is done.' }),
    innerEdgeFactor: PD.Numeric(1.5, { min: 0, max: 3, step: 0.1 }, { description: 'Factor to multiply the inner edge color with - for added contrast.' }),
};
export type MarkingProps = PD.Values<typeof MarkingParams>

export const SingleSample: number[][] = [[0, 0]];

/** Whether marked objects are tinted in the material color instead of by the marking pass. */
export function isMaterialColorMarker(rendererProps: Pick<RendererProps, 'colorMarker' | 'dimStrength'>, markingProps: MarkingProps) {
    // dimming unmarked objects is only supported in the material color
    return rendererProps.colorMarker && (!markingProps.enabled || rendererProps.dimStrength > 0);
}

type Props = {
    marking: MarkingProps;
    postprocessing: PostprocessingProps;
}

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

    /** jitter offsets the layer is accumulated from */
    private offsets = SingleSample;
    /** number of offsets the layer holds samples for, 0 if it is not valid */
    private sampleCount = 0;

    /** image of the last render to the drawing buffer, before marking was blended over it */
    private base: RenderTarget | undefined = undefined;

    constructor(private webgl: WebGLContext, width: number, height: number) {
        const { colorBufferFloat, textureFloat, colorBufferHalfFloat, textureHalfFloat } = webgl.extensions;
        this.depthTarget = webgl.createRenderTarget(width, height);
        // linear so that it can be used as antialiasing input
        this.maskTarget = webgl.createRenderTarget(width, height, true, 'uint8', 'linear');
        this.edgesTarget = webgl.createRenderTarget(width, height);
        const layerType = colorBufferHalfFloat && textureHalfFloat ? 'fp16' :
            colorBufferFloat && textureFloat ? 'float32' : 'uint8';
        this.layerTarget = webgl.createRenderTarget(width, height, false, layerType);

        this.edge = getEdgeRenderable(webgl, this.maskTarget.texture);
        this.overlay = getOverlayRenderable(webgl, this.edgesTarget.texture, this.maskTarget.texture);
        this.composite = createCopyRenderable(webgl, this.layerTarget.texture);
    }

    getByteCount() {
        return this.depthTarget.getByteCount() + this.maskTarget.getByteCount() + this.edgesTarget.getByteCount() + this.layerTarget.getByteCount() + (this.antialiasing?.getByteCount() ?? 0);
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
            this.sampleCount = 0;
            this.base = undefined;

            ValueCell.update(this.edge.values.uTexSizeInv, Vec2.set(this.edge.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.overlay.values.uTexSizeInv, Vec2.set(this.overlay.values.uTexSizeInv.ref.value, 1 / width, 1 / height));
            ValueCell.update(this.composite.values.uTexSize, Vec2.set(this.composite.values.uTexSize.ref.value, width, height));
        }
    }

    private update(props: MarkingProps, rendererProps: RendererProps) {
        const { highlightEdgeColor, selectEdgeColor, edgeScale, innerEdgeFactor, ghostEdgeStrength, highlightEdgeStrength, selectEdgeStrength } = props;

        const { values: edgeValues } = this.edge;
        ValueCell.updateIfChanged(edgeValues.uEdgeScale, Math.max(1, edgeScale * this.webgl.pixelRatio));

        const { values: overlayValues } = this.overlay;
        ValueCell.update(overlayValues.uHighlightEdgeColor, Color.toVec3Normalized(overlayValues.uHighlightEdgeColor.ref.value, highlightEdgeColor));
        ValueCell.update(overlayValues.uSelectEdgeColor, Color.toVec3Normalized(overlayValues.uSelectEdgeColor.ref.value, selectEdgeColor));
        ValueCell.updateIfChanged(overlayValues.uInnerEdgeFactor, innerEdgeFactor);
        ValueCell.updateIfChanged(overlayValues.uGhostEdgeStrength, ghostEdgeStrength);
        ValueCell.updateIfChanged(overlayValues.uHighlightEdgeStrength, highlightEdgeStrength);
        ValueCell.updateIfChanged(overlayValues.uSelectEdgeStrength, selectEdgeStrength);

        const { colorMarker, highlightColor, selectColor, highlightStrength, selectStrength } = rendererProps;
        const fill = colorMarker && !isMaterialColorMarker(rendererProps, props);
        ValueCell.update(overlayValues.uHighlightFillColor, Color.toVec3Normalized(overlayValues.uHighlightFillColor.ref.value, highlightColor));
        ValueCell.update(overlayValues.uSelectFillColor, Color.toVec3Normalized(overlayValues.uSelectFillColor.ref.value, selectColor));
        ValueCell.updateIfChanged(overlayValues.uHighlightFillStrength, fill ? highlightStrength : 0);
        ValueCell.updateIfChanged(overlayValues.uSelectFillStrength, fill ? selectStrength : 0);
    }

    /** Forgets the base image, needed when the last render did not go through `present`. */
    invalidate() {
        this.base = undefined;
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
    present(ctx: RenderContext, props: Props, base: RenderTarget, toDrawingBuffer: boolean, offsets: number[][], restart: boolean, samples: number) {
        const { viewport } = ctx.camera;
        const hasMarking = this.updateLayer(ctx, props, offsets, restart, samples);

        if (toDrawingBuffer) this.copyToDrawingBuffer(base, viewport);
        else if (hasMarking) base.bind();
        if (hasMarking) this.compositeLayer(viewport);

        this.base = toDrawingBuffer && !(ctx.camera instanceof StereoCamera) ? base : undefined;
    }

    /**
     * Redraws marking over the image of the last render to the drawing buffer, without
     * re-rendering the scene. Returns false if that image is not available.
     */
    redraw(ctx: RenderContext, props: Props, offsets: number[][], restart: boolean, samples: number): boolean {
        if (!this.base) return false;

        if (isTimingMode) this.webgl.timer.mark('MarkingPass.redraw');
        this.present(ctx, props, this.base, true, offsets, restart, samples);
        this.webgl.gl.flush();
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.redraw');
        return true;
    }

    /** Returns false if there is no marking. */
    private updateLayer(ctx: RenderContext, props: Props, offsets: number[][], restart: boolean, samples: number): boolean {
        const { renderer, camera, scene, frame } = ctx;
        if (camera instanceof StereoCamera || camera.disabled || !MarkingPass.hasMarking(scene, props)) {
            this.sampleCount = 0;
            return false;
        }
        if (restart || offsets !== this.offsets) {
            this.offsets = offsets;
            this.sampleCount = 0;
        }
        const end = Math.min(this.sampleCount + samples, offsets.length);
        if (this.sampleCount >= end) return true;

        if (isTimingMode) this.webgl.timer.mark('MarkingPass.updateLayer');
        const jitter = offsets.length > 1 && camera instanceof Camera;
        const { x, y, width, height } = camera.viewport;
        const depthTest = props.marking.ghostEdgeStrength < 1;

        renderer.setDrawingBufferSize(this.maskTarget.getWidth(), this.maskTarget.getHeight());
        renderer.setPixelRatio(this.webgl.pixelRatio);
        renderer.setViewport(x, y, width, height);
        this.update(props.marking, renderer.props);

        for (; this.sampleCount < end; ++this.sampleCount) {
            if (jitter) setJitter(camera, offsets[this.sampleCount]);
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

            this.renderSample(camera, props.postprocessing, 1 / (this.sampleCount + 1));
        }

        if (jitter) clearJitter(camera);
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.updateLayer');
        return true;
    }

    private renderMaskAa(camera: ICamera, props: PostprocessingProps): Texture {
        // no sharpening, it would break decoding the coverage from the mask
        if (!props.enabled || props.antialiasing.name === 'off') return this.maskTarget.texture;

        if (!this.antialiasing) {
            // linear, as the edge pass samples the mask between texels
            this.antialiasing = new AntialiasingPass(this.webgl, this.maskTarget.getWidth(), this.maskTarget.getHeight(), 'linear');
        }
        const { target } = this.antialiasing;
        return this.antialiasing.renderAntialiasingOnly(camera, this.maskTarget.texture, target, props)
            ? target.texture
            : this.maskTarget.texture;
    }

    /** Blends the marking of the current mask into the layer, with `weight` for the new sample. */
    private renderSample(camera: ICamera, postprocessingProps: PostprocessingProps, weight: number) {
        if (isTimingMode) this.webgl.timer.mark('MarkingPass.renderSample');
        const { gl, state } = this.webgl;
        const maskTexture = this.renderMaskAa(camera, postprocessingProps);

        if (this.edge.values.tMaskTexture.ref.value !== maskTexture) {
            ValueCell.update(this.edge.values.tMaskTexture, maskTexture);
            this.edge.update();
            ValueCell.update(this.overlay.values.tMaskTexture, maskTexture);
            this.overlay.update();
        }

        const viewport = camera.viewport;
        this.edgesTarget.bind();
        this.setEdgeState(viewport);
        this.edge.render();

        this.layerTarget.bind();
        this.setViewport(viewport);
        state.enable(gl.BLEND);
        // running average: the first sample (weight 1) replaces the previous content
        state.blendColor(0, 0, 0, weight);
        state.blendFunc(gl.CONSTANT_ALPHA, gl.ONE_MINUS_CONSTANT_ALPHA);
        state.blendEquation(gl.FUNC_ADD);
        this.overlay.render();
        if (isTimingMode) this.webgl.timer.markEnd('MarkingPass.renderSample');
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

const OverlaySchema = {
    ...QuadSchema,
    tEdgeTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    tMaskTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uTexSizeInv: UniformSpec('v2'),
    uHighlightEdgeColor: UniformSpec('v3'),
    uSelectEdgeColor: UniformSpec('v3'),
    uHighlightEdgeStrength: UniformSpec('f'),
    uSelectEdgeStrength: UniformSpec('f'),
    uGhostEdgeStrength: UniformSpec('f'),
    uInnerEdgeFactor: UniformSpec('f'),
    uHighlightFillColor: UniformSpec('v3'),
    uSelectFillColor: UniformSpec('v3'),
    uHighlightFillStrength: UniformSpec('f'),
    uSelectFillStrength: UniformSpec('f'),
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
        uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
        uHighlightEdgeColor: ValueCell.create(Vec3()),
        uSelectEdgeColor: ValueCell.create(Vec3()),
        uHighlightEdgeStrength: ValueCell.create(1),
        uSelectEdgeStrength: ValueCell.create(1),
        uGhostEdgeStrength: ValueCell.create(0),
        uInnerEdgeFactor: ValueCell.create(0),
        uHighlightFillColor: ValueCell.create(Vec3()),
        uSelectFillColor: ValueCell.create(Vec3()),
        uHighlightFillStrength: ValueCell.create(0),
        uSelectFillStrength: ValueCell.create(0),
    };

    const schema = { ...OverlaySchema };
    const renderItem = createComputeRenderItem(ctx, 'triangles', OverlayShaderCode, schema, values);

    return createComputeRenderable(renderItem, values);
}