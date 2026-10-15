/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import {
  type CopyRenderable,
  createCopyRenderable,
  getSharedCopyRenderable,
  QuadSchema,
  QuadValues,
} from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { DefineSpec, TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Viewport } from '../camera/util.js';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import { Color } from '@molstar/core/util/color';
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
  ghostEdgeStrength: PD.Numeric(
    0.3,
    { min: 0, max: 1, step: 0.1 },
    {
      description:
        'Opacity of the hidden edges that are covered by other geometry. When set to 1 and not dimming, one less geometry render pass is done, but hidden parts are then also tinted.',
    },
  ),
  innerEdgeFactor: PD.Numeric(
    1.5,
    { min: 0, max: 3, step: 0.1 },
    { description: 'Factor to multiply the inner edge color with - for added contrast.' },
  ),
};
export type MarkingProps = PD.Values<typeof MarkingParams>;

export const SingleSample: number[][] = [[0, 0]];

/** Whether marked objects are tinted and unmarked ones dimmed in the material color instead of by the marking pass. */
export function isMaterialColorMarker(rendererProps: Pick<RendererProps, 'colorMarker'>, markingProps: MarkingProps) {
  return rendererProps.colorMarker && !markingProps.enabled;
}

type Props = {
  marking: MarkingProps;
  postprocessing: PostprocessingProps;
  transparentBackground: boolean;
};

export type MarkingPresentOptions = {
  /** image to blend the marking over */
  base: RenderTarget;
  toDrawingBuffer: boolean;
  /** jitter offsets the layer is accumulated from */
  offsets: number[][];
  /** start the layer over instead of adding to it */
  restart: boolean;
  /** maximum number of samples to add */
  samples: number;
  /** how the base image was shaded, kept on tinted and dimmed objects */
  shading: MarkingShading | null;
};

export type MarkingShading =
  | { name: 'ssao'; ssao: Texture }
  /** the (denoised) base image relative to the direct `shaded` color gives the path traced occlusion and shadows */
  | { name: 'traced'; shaded: Texture };

type RenderContext = BaseRenderContext<ICamera | StereoCamera>;

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
  private dim: { target: RenderTarget; aaTarget: RenderTarget; renderable: DimRenderable } | null = null;

  /** jitter offsets the layer is accumulated from */
  private offsets = SingleSample;
  /** number of offsets the layer holds samples for, 0 if it is not valid */
  private sampleCount = 0;

  /** image of the last render to the drawing buffer, before marking was blended over it */
  private base: RenderTarget | undefined = undefined;
  /** shading of `base` */
  private baseShading: MarkingShading | null = null;

  constructor(
    private graphics: GraphicsContext,
    width: number,
    height: number,
  ) {
    const { floatRenderTarget, halfFloatRenderTarget } = graphics.capabilities;
    this.depthTarget = graphics.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'nearest');
    // linear so that it can be used as antialiasing input
    this.maskTarget = graphics.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'linear');
    this.edgesTarget = graphics.createRenderTarget(width, height);
    const layerType = halfFloatRenderTarget ? 'fp16' : floatRenderTarget ? 'float32' : 'uint8';
    this.layerTarget = graphics.createRenderTarget(width, height, 'none', layerType);

    this.edge = getEdgeRenderable(graphics, this.maskTarget.texture);
    this.overlay = getOverlayRenderable(graphics, this.edgesTarget.texture, this.maskTarget.texture);
    this.composite = createCopyRenderable(graphics, this.layerTarget.texture);
  }

  getByteCount() {
    const dimByteCount = this.dim ? this.dim.target.getByteCount() + this.dim.aaTarget.getByteCount() : 0;
    return (
      this.depthTarget.getByteCount() +
      this.maskTarget.getByteCount() +
      this.edgesTarget.getByteCount() +
      this.layerTarget.getByteCount() +
      (this.antialiasing?.getByteCount() ?? 0) +
      dimByteCount
    );
  }

  private setEdgeState(viewport: Viewport) {
    this.graphics.setRenderState(RenderState.Additive);
    this.setViewport(viewport);
    this.graphics.clear({ color: [0, 0, 0, 0] });
  }

  private setViewport(viewport: Viewport) {
    const { x, y, width, height } = viewport;
    this.graphics.setViewport(x, y, width, height);
    this.graphics.setScissor(x, y, width, height);
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
        ValueCell.update(
          this.dim.renderable.values.uTexSizeInv,
          Vec2.set(this.dim.renderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height),
        );
      }
      this.sampleCount = 0;
      this.invalidate();

      ValueCell.update(
        this.edge.values.uTexSizeInv,
        Vec2.set(this.edge.values.uTexSizeInv.ref.value, 1 / width, 1 / height),
      );
      ValueCell.update(
        this.overlay.values.uTexSizeInv,
        Vec2.set(this.overlay.values.uTexSizeInv.ref.value, 1 / width, 1 / height),
      );
      ValueCell.update(
        this.composite.values.uTexSize,
        Vec2.set(this.composite.values.uTexSize.ref.value, width, height),
      );
    }
  }

  private update(
    props: Props,
    rendererProps: RendererProps,
    dim: boolean,
    depthTest: boolean,
    shading: MarkingShading | null,
    base: Texture,
  ) {
    const {
      highlightEdgeColor,
      selectEdgeColor,
      edgeScale,
      innerEdgeFactor,
      ghostEdgeStrength,
      highlightEdgeStrength,
      selectEdgeStrength,
    } = props.marking;

    const { values: edgeValues } = this.edge;
    ValueCell.updateIfChanged(edgeValues.uEdgeScale, Math.max(1, edgeScale * this.graphics.pixelRatio));

    const { values: overlayValues } = this.overlay;
    ValueCell.update(
      overlayValues.uHighlightEdgeColor,
      Color.toVec3Normalized(overlayValues.uHighlightEdgeColor.ref.value, highlightEdgeColor),
    );
    ValueCell.update(
      overlayValues.uSelectEdgeColor,
      Color.toVec3Normalized(overlayValues.uSelectEdgeColor.ref.value, selectEdgeColor),
    );
    ValueCell.updateIfChanged(overlayValues.uInnerEdgeFactor, innerEdgeFactor);
    ValueCell.updateIfChanged(overlayValues.uGhostEdgeStrength, ghostEdgeStrength);
    ValueCell.updateIfChanged(overlayValues.uDepthTest, depthTest);
    ValueCell.updateIfChanged(overlayValues.uHighlightEdgeStrength, highlightEdgeStrength);
    ValueCell.updateIfChanged(overlayValues.uSelectEdgeStrength, selectEdgeStrength);

    const {
      colorMarker,
      highlightColor,
      selectColor,
      highlightStrength,
      selectStrength,
      dimColor,
      dimStrength,
      backgroundColor,
    } = rendererProps;
    const fill = colorMarker && !isMaterialColorMarker(rendererProps, props.marking);
    ValueCell.update(
      overlayValues.uHighlightFillColor,
      Color.toVec3Normalized(overlayValues.uHighlightFillColor.ref.value, highlightColor),
    );
    ValueCell.update(
      overlayValues.uSelectFillColor,
      Color.toVec3Normalized(overlayValues.uSelectFillColor.ref.value, selectColor),
    );
    ValueCell.updateIfChanged(overlayValues.uHighlightFillStrength, fill ? highlightStrength : 0);
    ValueCell.updateIfChanged(overlayValues.uSelectFillStrength, fill ? selectStrength : 0);
    ValueCell.update(overlayValues.uDimColor, Color.toVec3Normalized(overlayValues.uDimColor.ref.value, dimColor));
    ValueCell.updateIfChanged(overlayValues.uDimStrength, dim ? dimStrength : 0);

    const occlusionProps = props.postprocessing.occlusion;
    const shaded = dim || (fill && (highlightStrength > 0 || selectStrength > 0));
    const markingShading =
      !shaded || !shading
        ? 'off'
        : shading.name === 'traced'
          ? 'traced'
          : occlusionProps.name === 'on'
            ? 'ssao'
            : 'off';
    let needsUpdate = false;
    if (overlayValues.dMarkingShading.ref.value !== markingShading) {
      ValueCell.update(overlayValues.dMarkingShading, markingShading);
      needsUpdate = true;
    }
    if (shading?.name === 'ssao' && overlayValues.tSsaoDepth.ref.value !== shading.ssao) {
      ValueCell.update(overlayValues.tSsaoDepth, shading.ssao);
      needsUpdate = true;
    }
    if (
      shading?.name === 'traced' &&
      (overlayValues.tBase.ref.value !== base || overlayValues.tShaded.ref.value !== shading.shaded)
    ) {
      ValueCell.update(overlayValues.tBase, base);
      ValueCell.update(overlayValues.tShaded, shading.shaded);
      needsUpdate = true;
    }
    if (needsUpdate) this.overlay.update();

    if (markingShading === 'ssao' && occlusionProps.name === 'on') {
      ValueCell.update(
        overlayValues.uOcclusionColor,
        Color.toVec3Normalized(overlayValues.uOcclusionColor.ref.value, occlusionProps.params.color),
      );
    }
    if (markingShading !== 'off') {
      ValueCell.update(
        overlayValues.uFogColor,
        Color.toVec3Normalized(overlayValues.uFogColor.ref.value, backgroundColor),
      );
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
    else if (hasMarking) this.graphics.beginRenderPass(base, { label: 'marking composite' });
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

    if (isTimingMode) this.graphics.timer.group('MarkingPass.redraw');
    this.present(ctx, props, {
      base: this.base,
      toDrawingBuffer: true,
      offsets,
      restart,
      samples,
      shading: this.baseShading,
    });
    this.graphics.flush();
    if (isTimingMode) this.graphics.timer.groupEnd('MarkingPass.redraw');
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

    if (isTimingMode) this.graphics.timer.group('MarkingPass.updateLayer');
    const jitter = offsets.length > 1 && camera instanceof Camera;
    const { x, y, width, height } = camera.viewport;

    renderer.setDrawingBufferSize(this.maskTarget.getWidth(), this.maskTarget.getHeight());
    renderer.setPixelRatio(this.graphics.pixelRatio);
    renderer.setViewport(x, y, width, height);
    this.update(props, renderer.props, dim, depthTest, shading, base.texture);

    for (; this.sampleCount < end; ++this.sampleCount) {
      const offset = offsets[this.sampleCount];
      if (jitter) setJitter(camera, offset);
      renderer.update(camera, scene, frame);

      if (depthTest) {
        this.graphics.beginRenderPass(this.depthTarget, { label: 'marking depth' });
        renderer.clear(false, true);
        // when everything is marked there is nothing to occlude with, but the cleared depth is still needed
        if (scene.markerAverage !== 1) renderer.renderMarkingDepth(scene.primitives, camera);
      }

      this.graphics.beginRenderPass(this.maskTarget, { label: 'marking mask' });
      renderer.clear(false, true);
      renderer.renderMarkingMask(scene.primitives, camera, depthTest ? this.depthTarget.texture : null);

      // the scene SSAO is computed without jitter, see `MultiSamplePass`
      ValueCell.update(
        this.overlay.values.uOcclusionOffset,
        Vec2.set(
          this.overlay.values.uOcclusionOffset.ref.value,
          jitter ? offset[0] / width : 0,
          jitter ? offset[1] / height : 0,
        ),
      );
      this.renderSample(camera, props.postprocessing, 1 / (this.sampleCount + 1), dim);
    }

    if (jitter) clearJitter(camera);
    if (isTimingMode) this.graphics.timer.groupEnd('MarkingPass.updateLayer');
    return true;
  }

  private getAntialiasing(props: PostprocessingProps): AntialiasingPass | null {
    // no sharpening, it would break decoding the coverage from the mask
    if (!props.enabled || props.antialiasing.name === 'off') return null;

    if (!this.antialiasing) {
      // linear, as the edge pass samples the mask between texels
      this.antialiasing = new AntialiasingPass(
        this.graphics,
        this.maskTarget.getWidth(),
        this.maskTarget.getHeight(),
        'linear',
      );
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
        target: this.graphics.createRenderTarget(width, height, 'none', 'uint8', 'linear'),
        aaTarget: this.graphics.createRenderTarget(width, height, 'none'),
        renderable: getDimRenderable(this.graphics, this.depthTarget.texture, this.maskTarget.texture),
      };
    }
    return this.dim;
  }

  /** Renders the antialiased coverage of the visible unmarked geometry, from the depth and raw mask of the current sample. */
  private renderDim(camera: ICamera, postprocessingProps: PostprocessingProps): Texture {
    const { target, aaTarget, renderable } = this.getDim();
    const { values } = renderable;
    ValueCell.updateIfChanged(values.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
    ValueCell.updateIfChanged(values.uNear, camera.near);
    ValueCell.updateIfChanged(values.uFar, camera.far);
    ValueCell.updateIfChanged(values.uFogNear, camera.fogNear);
    ValueCell.updateIfChanged(values.uFogFar, camera.fogFar);

    this.graphics.beginRenderPass(target, { label: 'marking dim' });
    this.setViewport(camera.viewport);
    this.graphics.setRenderState(RenderState.Default);
    this.graphics.clear({ color: [0, 0, 0, 0] });
    renderable.render();

    return this.antialias(camera, postprocessingProps, target, aaTarget);
  }

  /** Blends the marking of the current mask into the layer, with `weight` for the new sample. */
  private renderSample(camera: ICamera, postprocessingProps: PostprocessingProps, weight: number, dim: boolean) {
    if (isTimingMode) this.graphics.timer.group('MarkingPass.renderSample');

    const dimTexture = dim ? this.renderDim(camera, postprocessingProps) : this.maskTarget.texture;
    const aa = this.getAntialiasing(postprocessingProps);
    const maskTexture = aa
      ? this.antialias(camera, postprocessingProps, this.maskTarget, aa.target)
      : this.maskTarget.texture;

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
    this.graphics.beginRenderPass(this.edgesTarget, { label: 'marking edges' });
    this.setEdgeState(viewport);
    this.edge.render();

    this.renderOverlay(viewport, weight);
    if (isTimingMode) this.graphics.timer.groupEnd('MarkingPass.renderSample');
  }

  /** Blends the overlay of the current sample's textures into the layer, with `weight` for the new sample. */
  private renderOverlay(viewport: Viewport, weight: number) {
    this.graphics.beginRenderPass(this.layerTarget, { label: 'marking overlay' });
    this.setViewport(viewport);
    // running average: the first sample (weight 1) replaces the previous content
    this.graphics.setBlendConstant(weight, weight, weight, weight);
    this.graphics.setRenderState(RenderState.Constant);
    this.overlay.render();
  }

  /** Blends the layer over the currently bound framebuffer. */
  private compositeLayer(viewport: Viewport) {
    this.graphics.setRenderState(RenderState.Premultiplied);
    this.setViewport(viewport);
    this.composite.render();
  }

  private copyToDrawingBuffer(src: RenderTarget, viewport: Viewport) {
    const copy = getSharedCopyRenderable(this.graphics, src.texture);
    this.graphics.beginRenderPass('drawing-buffer', { label: 'copy' });
    this.graphics.setRenderState(RenderState.Default);
    this.setViewport(viewport);
    copy.render();
  }
}

//

export const EdgeSchema = {
  ...QuadSchema,
  tMaskTexture: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  uTexSizeInv: UniformSpec('v2'),
  uEdgeScale: UniformSpec('f'),
};
export const EdgeShaderCode = ShaderCode('edge');
type EdgeRenderable = ComputeRenderable<Values<typeof EdgeSchema>>;

function getEdgeRenderable(ctx: GraphicsContext, maskTexture: Texture): EdgeRenderable {
  const width = maskTexture.getWidth();
  const height = maskTexture.getHeight();

  const values: Values<typeof EdgeSchema> = {
    ...QuadValues,
    tMaskTexture: ValueCell.create(maskTexture),
    uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),
    uEdgeScale: ValueCell.create(1),
  };

  const schema = { ...EdgeSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', EdgeShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

export const DimSchema = {
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
export const DimShaderCode = ShaderCode('dim');
type DimRenderable = ComputeRenderable<Values<typeof DimSchema>>;

function getDimRenderable(ctx: GraphicsContext, depthTexture: Texture, maskTexture: Texture): DimRenderable {
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
  const renderItem = ctx.createComputeRenderItem('triangles', DimShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

export const OverlaySchema = {
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
export const OverlayShaderCode = ShaderCode('overlay');
type OverlayRenderable = ComputeRenderable<Values<typeof OverlaySchema>>;

function getOverlayRenderable(ctx: GraphicsContext, edgeTexture: Texture, maskTexture: Texture): OverlayRenderable {
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
  const renderItem = ctx.createComputeRenderItem('triangles', OverlayShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
