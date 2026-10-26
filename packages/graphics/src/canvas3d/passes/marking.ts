/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { type PassBuilder } from './pass-builder.js';
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
import type { Renderbuffer } from '@molstar/graphics/gl/renderbuffer';
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
import { updateOcclusionUpsample } from './ssao.js';
import type { PassTargets } from './pass-targets.js';

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
  /** the renderer's props where they are not the live renderer's (preparing props a canvas is not set to) */
  renderer?: Partial<RendererProps>;
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
  /** draw as if all objects were marked, to prepare the pipelines of marking ahead (see `DrawPass.plan`) */
  asIfMarked?: boolean;
};

export type MarkingShading =
  /** the SSAO targets, upsampled depth-aware with the scene's opaque `depth` (packed without depth textures) */
  | { name: 'ssao'; ssao: Texture; depth: Texture; packedDepth: boolean }
  /** the (denoised) base image relative to the direct `shaded` color gives the path traced occlusion and shadows */
  | { name: 'traced'; shaded: Texture };

type RenderContext = BaseRenderContext<ICamera | StereoCamera>;

export class MarkingPass {
  private readonly graphics: GraphicsContext;
  static isEnabled(props: MarkingProps) {
    return props.enabled;
  }

  static hasMarking(scene: Scene, props: { marking: MarkingProps }) {
    return MarkingPass.isEnabled(props.marking) && scene.markerAverage > 0;
  }

  private readonly depthTarget: RenderTarget;
  private readonly maskTarget: RenderTarget;
  /** premultiplied marking, the running average of the samples taken so far */
  private readonly layerTarget: RenderTarget;

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

  private readonly renderBuilder: PassBuilder;

  /**
   * `depthStencil`: the depth-stencil buffer of the passes' size the depth and the mask are drawn with (depth cleared
   * per pass, stencil per solid interior cap), not read after; shared with passes that use it earlier in the frame
   */
  constructor(
    private readonly targets: PassTargets,
    width: number,
    height: number,
    depthStencil: Renderbuffer,
  ) {
    const graphics = (this.graphics = targets.graphics);
    const { floatBlend, halfFloatRenderTarget } = graphics.capabilities;
    this.depthTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'nearest');
    // linear so that it can be used as antialiasing input
    this.maskTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');
    depthStencil.attachFramebuffer(this.depthTarget.framebuffer);
    depthStencil.attachFramebuffer(this.maskTarget.framebuffer);
    const layerType = halfFloatRenderTarget ? 'fp16' : floatBlend ? 'float32' : 'uint8';
    this.layerTarget = targets.createRenderTarget(width, height, 'none', layerType);

    this.overlay = getOverlayRenderable(graphics, this.maskTarget.texture);
    this.composite = createCopyRenderable(graphics, this.layerTarget.texture);
    this.renderBuilder = targets.createRenderBuilder();
  }

  getByteCount() {
    const dimByteCount = this.dim ? this.dim.target.getByteCount() + this.dim.aaTarget.getByteCount() : 0;
    return (
      this.depthTarget.getByteCount() +
      this.maskTarget.getByteCount() +
      this.layerTarget.getByteCount() +
      (this.antialiasing?.getByteCount() ?? 0) +
      dimByteCount
    );
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
    camera: ICamera,
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

    const { values: overlayValues } = this.overlay;
    ValueCell.updateIfChanged(overlayValues.uEdgeScale, Math.max(1, edgeScale * this.graphics.pixelRatio));
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
    if (shading?.name === 'ssao') {
      if (overlayValues.tSsaoDepth.ref.value !== shading.ssao || overlayValues.tDepth.ref.value !== shading.depth) {
        ValueCell.update(overlayValues.tSsaoDepth, shading.ssao);
        ValueCell.update(overlayValues.tDepth, shading.depth);
        needsUpdate = true;
      }
      if (overlayValues.dPackedDepth.ref.value !== shading.packedDepth) {
        ValueCell.update(overlayValues.dPackedDepth, shading.packedDepth);
        needsUpdate = true;
      }
      const width = 1 / overlayValues.uTexSizeInv.ref.value[0];
      const height = 1 / overlayValues.uTexSizeInv.ref.value[1];
      if (updateOcclusionUpsample(overlayValues, shading.ssao, Math.round(width), Math.round(height)))
        needsUpdate = true;
      ValueCell.updateIfChanged(overlayValues.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
      ValueCell.updateIfChanged(overlayValues.uNear, camera.near);
      ValueCell.updateIfChanged(overlayValues.uFar, camera.far);
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
   * or in place. Preparing (`b.drawing` false) leaves the layer and the remembered base as they are.
   */
  present(b: PassBuilder, ctx: RenderContext, props: Props, options: MarkingPresentOptions) {
    const { base, toDrawingBuffer, shading } = options;
    const { viewport } = ctx.camera;
    const hasMarking = this.updateLayer(b, ctx, props, options);

    if (toDrawingBuffer) this.copyToDrawingBuffer(b, base, viewport);
    else if (hasMarking) b.begin(base, { label: 'marking composite' });
    if (hasMarking) this.compositeLayer(b, viewport);

    if (!b.drawing) return;
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
    this.present(this.renderBuilder, ctx, props, {
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
  private updateLayer(b: PassBuilder, ctx: RenderContext, props: Props, options: MarkingPresentOptions): boolean {
    const { base, offsets, restart, samples, shading, asIfMarked } = options;
    const { renderer, camera, scene, frame } = ctx;
    if (camera instanceof StereoCamera || camera.disabled || !(asIfMarked || MarkingPass.hasMarking(scene, props))) {
      if (b.drawing) this.sampleCount = 0;
      return false;
    }

    const colorMarker = props.renderer?.colorMarker ?? renderer.props.colorMarker;
    const dimStrength = props.renderer?.dimStrength ?? renderer.props.dimStrength;
    // nothing to dim when everything is marked
    const dim = colorMarker && dimStrength > 0 && scene.markerAverage < 1;
    // dimming needs to know where marked objects are visible
    const depthTest = props.marking.ghostEdgeStrength < 1 || dim;

    if (!b.drawing) {
      // the samples draw alike, one is enough to request their pipelines
      this.update(props, renderer.props, dim, depthTest, shading, base.texture, camera);
      this.planSample(b, ctx, props, depthTest, dim, !!asIfMarked, 1);
      return true;
    }

    if (restart || offsets !== this.offsets) {
      this.offsets = offsets;
      this.sampleCount = 0;
    }

    const end = Math.min(this.sampleCount + samples, offsets.length);
    if (this.sampleCount >= end) {
      // traced shading comes from the base image, which keeps converging; the textures of a single sample are still there
      if (shading?.name === 'traced' && offsets.length === 1) {
        this.update(props, renderer.props, dim, depthTest, shading, base.texture, camera);
        this.renderOverlay(b, ctx.camera.viewport, 1);
      }
      return true;
    }

    if (isTimingMode) this.graphics.timer.group('MarkingPass.updateLayer');
    const jitter = offsets.length > 1 && camera instanceof Camera;
    const { x, y, width, height } = camera.viewport;

    renderer.setDrawingBufferSize(this.maskTarget.getWidth(), this.maskTarget.getHeight());
    renderer.setPixelRatio(this.graphics.pixelRatio);
    renderer.setViewport(x, y, width, height);
    this.update(props, renderer.props, dim, depthTest, shading, base.texture, camera);

    for (; this.sampleCount < end; ++this.sampleCount) {
      const offset = offsets[this.sampleCount];
      if (jitter) setJitter(camera, offset);
      renderer.update(camera, scene, frame);

      // the scene SSAO is computed without jitter, see `MultiSamplePass`
      ValueCell.update(
        this.overlay.values.uOcclusionOffset,
        Vec2.set(
          this.overlay.values.uOcclusionOffset.ref.value,
          jitter ? offset[0] / width : 0,
          jitter ? -offset[1] / height : 0,
        ),
      );
      this.planSample(b, ctx, props, depthTest, dim, !!asIfMarked, 1 / (this.sampleCount + 1));
    }

    if (jitter) clearJitter(camera);
    if (isTimingMode) this.graphics.timer.groupEnd('MarkingPass.updateLayer');
    return true;
  }

  /** The depth and mask of a sample and its marking blended into the layer (`asIfMarked`: all objects) */
  private planSample(
    b: PassBuilder,
    ctx: RenderContext,
    props: Props,
    depthTest: boolean,
    dim: boolean,
    asIfMarked: boolean,
    weight: number,
  ) {
    const { renderer, scene } = ctx;
    const camera = ctx.camera as ICamera;
    // direct volumes write the depth where they look solid (their `surfaceOpacity`); none: no pass (a replay of an empty
    // list still costs)
    const volumes = this.graphics.capabilities.fragDepth && scene.volumes.renderables.length > 0;
    if (depthTest) {
      b.begin(this.depthTarget, { label: 'marking depth' });
      if (b.drawing) renderer.clearDepth(true);
      // when everything is marked there is nothing to occlude with, but the cleared depth is still needed
      if (asIfMarked || scene.markerAverage !== 1) {
        b.renderer(
          { name: 'markingDepth' },
          scene.primitives,
          () => renderer.renderMarkingDepth(scene.primitives, camera),
          asIfMarked,
        );
        if (volumes) {
          b.renderer(
            { name: 'markingDepth' },
            scene.volumes,
            () => renderer.renderMarkingDepth(scene.volumes, camera),
            asIfMarked,
          );
        }
      }
    }

    b.begin(this.maskTarget, { label: 'marking mask' });
    if (b.drawing) renderer.clear(false, true);
    // the mask draws marked objects only: as if marked, all objects
    b.renderer(
      { name: 'markingMask' },
      scene.primitives,
      () => renderer.renderMarkingMask(scene.primitives, camera, depthTest ? this.depthTarget.texture : null),
      asIfMarked,
    );
    if (volumes) {
      b.renderer(
        { name: 'markingMask' },
        scene.volumes,
        () => renderer.renderMarkingMask(scene.volumes, camera, depthTest ? this.depthTarget.texture : null),
        asIfMarked,
      );
    }

    this.renderSample(b, camera, props.postprocessing, weight, dim);
  }

  private getAntialiasing(props: PostprocessingProps): AntialiasingPass | null {
    // no sharpening, it would break decoding the coverage from the mask
    if (!props.enabled || props.antialiasing.name === 'off') return null;

    if (!this.antialiasing) {
      // linear, as the overlay samples the mask between texels for the edges
      this.antialiasing = new AntialiasingPass(
        this.targets,
        this.maskTarget.getWidth(),
        this.maskTarget.getHeight(),
        'linear',
      );
    }
    return this.antialiasing;
  }

  /** Antialiases `input` into `output`, returns the texture to use. */
  private antialias(
    b: PassBuilder,
    camera: ICamera,
    props: PostprocessingProps,
    input: RenderTarget,
    output: RenderTarget,
  ): Texture {
    const aa = this.getAntialiasing(props);
    return aa?.renderAntialiasingOnly(b, camera, input.texture, output, props) ? output.texture : input.texture;
  }

  private getDim() {
    if (!this.dim) {
      const width = this.maskTarget.getWidth();
      const height = this.maskTarget.getHeight();
      this.dim = {
        // linear so that it can be used as antialiasing input
        target: this.targets.createRenderTarget(width, height, 'none', 'uint8', 'linear'),
        aaTarget: this.targets.createRenderTarget(width, height, 'none'),
        renderable: getDimRenderable(this.graphics, this.depthTarget.texture, this.maskTarget.texture),
      };
    }
    return this.dim;
  }

  /** Renders the antialiased coverage of the visible unmarked geometry, from the depth and raw mask of the current sample. */
  private renderDim(b: PassBuilder, camera: ICamera, postprocessingProps: PostprocessingProps): Texture {
    const { target, aaTarget, renderable } = this.getDim();
    const { values } = renderable;
    ValueCell.updateIfChanged(values.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
    ValueCell.updateIfChanged(values.uNear, camera.near);
    ValueCell.updateIfChanged(values.uFar, camera.far);
    ValueCell.updateIfChanged(values.uFogNear, camera.fogNear);
    ValueCell.updateIfChanged(values.uFogFar, camera.fogFar);

    b.begin(target, { label: 'marking dim' });
    if (b.drawing) {
      this.setViewport(camera.viewport);
      this.graphics.setRenderState(RenderState.Default);
      this.graphics.clear({ color: [0, 0, 0, 0] });
    }
    b.quad(renderable, RenderState.Default);

    return this.antialias(b, camera, postprocessingProps, target, aaTarget);
  }

  /** Blends the marking of the current mask into the layer, with `weight` for the new sample. */
  private renderSample(
    b: PassBuilder,
    camera: ICamera,
    postprocessingProps: PostprocessingProps,
    weight: number,
    dim: boolean,
  ) {
    const timing = isTimingMode && b.drawing;
    if (timing) this.graphics.timer.group('MarkingPass.renderSample');

    const dimTexture = dim ? this.renderDim(b, camera, postprocessingProps) : this.maskTarget.texture;
    const aa = this.getAntialiasing(postprocessingProps);
    const maskTexture = aa
      ? this.antialias(b, camera, postprocessingProps, this.maskTarget, aa.target)
      : this.maskTarget.texture;

    if (this.overlay.values.tMaskTexture.ref.value !== maskTexture) {
      ValueCell.update(this.overlay.values.tMaskTexture, maskTexture);
      this.overlay.update();
    }
    if (this.overlay.values.tDimTexture.ref.value !== dimTexture) {
      ValueCell.update(this.overlay.values.tDimTexture, dimTexture);
      this.overlay.update();
    }

    this.renderOverlay(b, camera.viewport, weight);
    if (timing) this.graphics.timer.groupEnd('MarkingPass.renderSample');
  }

  /** Blends the overlay of the current sample's textures into the layer, with `weight` for the new sample. */
  private renderOverlay(b: PassBuilder, viewport: Viewport, weight: number) {
    b.begin(this.layerTarget, { label: 'marking overlay' });
    if (b.drawing) {
      this.setViewport(viewport);
      // running average: the first sample (weight 1) replaces the previous content
      this.graphics.setBlendConstant(weight, weight, weight, weight);
    }
    b.quad(this.overlay, RenderState.Constant);
  }

  /** Blends the layer over the currently bound framebuffer. */
  private compositeLayer(b: PassBuilder, viewport: Viewport) {
    if (b.drawing) this.setViewport(viewport);
    b.quad(this.composite, RenderState.Premultiplied);
  }

  private copyToDrawingBuffer(b: PassBuilder, src: RenderTarget, viewport: Viewport) {
    const copy = getSharedCopyRenderable(this.graphics, src.texture);
    b.begin('drawing-buffer', { label: 'copy' });
    if (b.drawing) this.setViewport(viewport);
    b.quad(copy, RenderState.Default);
  }
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
  uEdgeScale: UniformSpec('f'),
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
  uOcclusionOffset: { ...UniformSpec('v2'), windowRect: 'offset' as const },
  /** the SSAO targets have a lower resolution, see `getSsaoUpsampled`; with the scene's opaque depth */
  dOcclusionUpsample: DefineSpec('boolean'),
  uOcclusionTexSize: UniformSpec('v2'),
  tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dPackedDepth: DefineSpec('boolean'),
  uIsOrtho: UniformSpec('f'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  tBase: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tShaded: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dMarkingShading: DefineSpec('string', ['off', 'ssao', 'traced']),
};
export const OverlayShaderCode = ShaderCode('overlay');
type OverlayRenderable = ComputeRenderable<Values<typeof OverlaySchema>>;

function getOverlayRenderable(ctx: GraphicsContext, maskTexture: Texture): OverlayRenderable {
  const width = maskTexture.getWidth();
  const height = maskTexture.getHeight();

  const values: Values<typeof OverlaySchema> = {
    ...QuadValues,
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
    uEdgeScale: ValueCell.create(1),
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
    dOcclusionUpsample: ValueCell.create(false),
    uOcclusionTexSize: ValueCell.create(Vec2()),
    tDepth: ValueCell.create(maskTexture),
    dPackedDepth: ValueCell.create(false),
    uIsOrtho: ValueCell.create(0),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    tBase: ValueCell.create(maskTexture),
    tShaded: ValueCell.create(maskTexture),
    dMarkingShading: ValueCell.create('off'),
  };

  const schema = { ...OverlaySchema };
  const renderItem = ctx.createComputeRenderItem('triangles', OverlayShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
