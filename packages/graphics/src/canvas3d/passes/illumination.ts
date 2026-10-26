/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { type PassBuilder } from './pass-builder.js';
import { type CopyRenderable, QuadSchema, QuadValues, createCopyRenderable } from '@molstar/graphics/gl/compute/util';
import { DefineSpec, TextureSpec, UniformSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import type { Texture } from '@molstar/graphics/gl/texture';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { ValueCell } from '@molstar/core/util';
import { isDebugMode, isTimingMode } from '@molstar/core/util/debug';
import type { Renderer, RendererProps } from '@molstar/graphics/gl/renderer';
import { Camera, type ICamera } from '../camera.js';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Color } from '@molstar/core/util/color/color';
import { AntialiasingPass, PostprocessingPass, type PostprocessingProps } from './postprocessing.js';
import type { DrawPass } from './draw.js';
import { type MarkingProps, type MarkingShading, SingleSample } from './marking.js';
import { DofPass } from './dof.js';
import { TraceShaderCode, TracingParams, TracingPass } from './tracing.js';
import { getJitterOffsets, setJitter, clearJitter } from './jitter.js';
import type { MultiSampleProps } from './multi-sample.js';
import { clamp, lerp } from '@molstar/core/math/interpolate';
import { type SsaoProps, updateOcclusionUpsample } from './ssao.js';
import { OutlinePass } from './outline.js';
import { BloomPass } from './bloom.js';
import type { RenderContext } from '../util.js';
import { Denoiser, isComputeDenoiseSupported } from '@molstar/graphics/gl/compute/denoise';

let IlluminationWarningShown = false;

const IlluminationCapabilities = ['floatBlend', 'depthTexture', 'multipleRenderTargets'] as const;

/** Why illumination is not available, if it is not: programs or capabilities missing */
function getIlluminationUnsupportedReason(graphics: GraphicsContext): string | undefined {
  if (!graphics.isProgramSupported(TraceShaderCode) || !graphics.isProgramSupported(ComposeShaderCode)) {
    return 'programs not supported';
  }
  const missing = IlluminationCapabilities.filter((k) => !graphics.capabilities[k]);
  if (missing.length) return `missing ${missing.join(', ')}`;
}

function checkIlluminationSupport(graphics: GraphicsContext) {
  const reason = getIlluminationUnsupportedReason(graphics);
  if (reason && isDebugMode && !IlluminationWarningShown) {
    console.log(`"illumination" not supported: ${reason}`);
    IlluminationWarningShown = true;
  }
  return !reason;
}

type Props = {
  transparentBackground: boolean;
  dpoitIterations: number;
  volumeResolutionScale: number;
  illumination: IlluminationProps;
  renderer: RendererProps;
  postprocessing: PostprocessingProps;
  marking: MarkingProps;
  multiSample: MultiSampleProps;
};

export const IlluminationParams = {
  enabled: PD.Boolean(false),
  maxIterations: PD.Numeric(
    5,
    { min: 0, max: 16, step: 1 },
    { description: 'Maximum number of tracing iterations. Final iteration count is 2^x.' },
  ),
  denoise: PD.Boolean(true),
  denoiseThreshold: PD.Interval(
    [0.15, 1],
    { min: 0, max: 4, step: 0.01 },
    {
      description: 'Threshold for denoising. Automatically adjusted within given interval based on current iteration.',
    },
  ),
  ignoreOutline: PD.Boolean(true, {
    description:
      'Ignore outline in illumination pass where it is generally not needed for visual clarity. Useful when illumination is often toggled on/off.',
  }),
  ...TracingParams,
};
export type IlluminationProps = PD.Values<typeof IlluminationParams>;

export class IlluminationPass {
  private readonly tracing: TracingPass;

  private readonly transparentTarget: RenderTarget;
  private readonly outputTarget: RenderTarget;

  readonly packedDepth: boolean;

  private readonly copyRenderable: CopyRenderable;
  private readonly composeRenderable: ComposeRenderable;

  /** the weighted average of the samples so far, see `renderMultiSample` */
  private multiSampleComposeTarget: RenderTarget;
  /** the first sample */
  private multiSampleHoldTarget: RenderTarget;
  private multiSampleCompose: MultiSampleComposeRenderable;
  private readonly renderBuilder: PassBuilder;

  private _iteration = 0;
  get iteration() {
    return this._iteration;
  }

  private _colorTarget: RenderTarget;
  get colorTarget() {
    return this._colorTarget;
  }

  private markingShading: MarkingShading;
  /** the denoiser as a kernel (compute shaders), created on first use; else the compose program denoises */
  private denoiser: Denoiser | undefined;

  private _supported = false;
  get supported() {
    return this._supported;
  }

  getByteCount() {
    if (!this._supported) return 0;
    return (
      this.tracing.getByteCount() +
      this.transparentTarget.getByteCount() +
      this.outputTarget.getByteCount() +
      this.multiSampleComposeTarget.getByteCount() +
      this.multiSampleHoldTarget.getByteCount()
    );
  }

  getMaxIterations(props: IlluminationProps) {
    return Math.pow(2, props.maxIterations);
  }

  static getUnsupportedReason(graphics: GraphicsContext) {
    return getIlluminationUnsupportedReason(graphics);
  }

  static isSupported(graphics: GraphicsContext) {
    return checkIlluminationSupport(graphics);
  }

  static isEnabled(graphics: GraphicsContext, props: IlluminationProps) {
    return props.enabled && checkIlluminationSupport(graphics);
  }

  constructor(
    private readonly graphics: GraphicsContext,
    private readonly drawPass: DrawPass,
  ) {
    if (!IlluminationPass.isSupported(graphics)) return;

    const { colorTarget } = drawPass;
    const width = colorTarget.getWidth();
    const height = colorTarget.getHeight();

    this.tracing = new TracingPass(graphics, this.drawPass);
    this.markingShading = { name: 'traced', shaded: this.tracing.shadedTextureOpaque };

    this.transparentTarget = drawPass.targets.createRenderTarget(width, height, 'none', 'uint8', 'nearest');
    this.outputTarget = drawPass.targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');

    this.copyRenderable = createCopyRenderable(graphics, this.transparentTarget.texture);

    this.composeRenderable = getComposeRenderable(
      graphics,
      this.tracing.accumulateTarget.texture,
      this.tracing.normalTextureOpaque,
      this.tracing.colorTextureOpaque,
      this.drawPass.depthTextureOpaque,
      this.drawPass.depthTargetTransparent.texture,
      this.drawPass.postprocessing.outline.target.texture,
      this.transparentTarget.texture,
      this.drawPass.postprocessing.ssao.ssaoDepthTexture,
      this.drawPass.postprocessing.ssao.ssaoDepthTransparentTexture,
      this.drawPass.postprocessing.bloom.compositeTarget.texture,
      false,
    );

    this.multiSampleComposeTarget = drawPass.targets.createRenderTarget(width, height, 'none', 'float32');
    this.multiSampleHoldTarget = drawPass.targets.createRenderTarget(width, height, 'none');
    this.multiSampleCompose = getMultiSampleComposeRenderable(graphics, this.outputTarget.texture);

    this.renderBuilder = drawPass.targets.createRenderBuilder();
    this._supported = true;
  }

  private renderInput(b: PassBuilder, renderer: Renderer, camera: ICamera, scene: Scene, props: Props) {
    const timing = isTimingMode && b.drawing;
    if (timing) this.graphics.timer.group('IlluminationPass.renderInput');

    const hasTransparent = scene.opacityAverage < 1 || scene.volumes.renderables.length > 0;
    const depthOpaque = this.drawPass.depthTextureOpaque;

    b.begin(this.transparentTarget, {
      label: 'transparent clear',
      clear: { color: [0, 0, 0, 0] },
    });

    const outlineEnabled =
      PostprocessingPass.isTransparentOutlineEnabled(props.postprocessing) && !props.illumination.ignoreOutline;
    const dofEnabled = DofPass.isEnabled(props.postprocessing);
    const ssaoEnabled = PostprocessingPass.isTransparentSsaoEnabled(scene, props.postprocessing);
    const bloomEnabled = BloomPass.isActive(scene, props.postprocessing);

    if ((outlineEnabled || dofEnabled || ssaoEnabled || bloomEnabled) && b.drawing) {
      b.begin(this.drawPass.depthTargetTransparent, { label: 'transparent depth clear' });
      renderer.clearDepth(true);
    }

    if (hasTransparent) {
      // volumes under the transparent primitives (one layer, not split at their surfaces)
      if (scene.volumes.renderables.length > 0) {
        b.begin(this.transparentTarget, { label: 'volume' });
        b.renderer({ name: 'volume' }, scene.volumes, () => renderer.renderVolume(scene.volumes, camera, depthOpaque));
      }

      if (this.drawPass.transparency === 'wboit') {
        this.drawPass.wboit.bind(b);
        b.renderer({ name: 'wboitTransparent' }, scene.primitives, () =>
          renderer.renderWboitTransparent(scene.primitives, camera, depthOpaque),
        );

        b.begin(this.transparentTarget, { label: 'wboit evaluate' });
        this.drawPass.wboit.render(b);
      } else if (this.drawPass.transparency === 'dpoit') {
        const dpoitTextures = this.drawPass.dpoit.bind(b);
        b.renderer({ name: 'dpoitTransparent' }, scene.primitives, () =>
          renderer.renderDpoitTransparent(scene.primitives, camera, depthOpaque, dpoitTextures),
        );

        // preparing: one peeling layer draws as all of them
        const iterations = props.dpoitIterations;
        const layers = b.drawing ? iterations : Math.min(iterations, 1);
        for (let i = 0; i < layers; i++) {
          if (timing) this.graphics.timer.group('DpoitPass.layer');
          const dpoitTextures = this.drawPass.dpoit.bindDualDepthPeeling(b);
          b.renderer({ name: 'dpoitTransparent' }, scene.primitives, () =>
            renderer.renderDpoitTransparent(scene.primitives, camera, depthOpaque, dpoitTextures),
          );

          if (iterations > 1) {
            b.begin(this.transparentTarget, { label: 'dpoit blend back' });
            this.drawPass.dpoit.renderBlendBack(b);
          }
          if (timing) this.graphics.timer.groupEnd('DpoitPass.layer');
        }

        // evaluate dpoit
        b.begin(this.transparentTarget, { label: 'dpoit evaluate' });
        this.drawPass.dpoit.render(b);
      } else {
        depthOpaque.attachFramebuffer(this.transparentTarget.framebuffer, 'depth');
        b.begin(this.transparentTarget, { label: 'transparent' });
        b.renderer({ name: 'blendedTransparent' }, scene.primitives, () =>
          renderer.renderBlendedTransparent(scene.primitives, camera),
        );
        depthOpaque.detachFramebuffer(this.transparentTarget.framebuffer, 'depth');
      }

      if (outlineEnabled || dofEnabled || ssaoEnabled || bloomEnabled) {
        b.begin(this.drawPass.depthTargetTransparent, { label: 'transparent depth' });
        if (scene.opacityAverage < 1) {
          b.renderer({ name: 'depthTransparent' }, scene.primitives, () =>
            renderer.renderDepthTransparent(scene.primitives, camera, depthOpaque),
          );
        }
      }

      if (ssaoEnabled) {
        this.drawPass.postprocessing.ssao.update(
          camera,
          scene,
          props.postprocessing.occlusion.params as SsaoProps,
          true,
        );
        this.drawPass.postprocessing.ssao.render(b, camera);
      }
    }

    //

    b.begin(this.tracing.composeTarget, {
      label: 'compose clear',
      clear: { color: [0, 0, 0, 0] },
    });
    if (timing) this.graphics.timer.groupEnd('IlluminationPass.renderInput');
  }

  private hasHelpers(helper: RenderContext<Camera>['helper']) {
    return helper.debug.isEnabled || helper.pointer.isEnabled || helper.handle.isEnabled || helper.camera.isEnabled;
  }

  private renderHelpers(b: PassBuilder, ctx: RenderContext<Camera>) {
    const { renderer, camera, helper, frame } = ctx;
    b.begin(this.transparentTarget, { label: 'helpers' });
    if (helper.debug.isEnabled || helper.pointer.isEnabled) {
      this.drawPass.depthTextureOpaque.attachFramebuffer(this.transparentTarget.framebuffer, 'depth');
      b.begin(this.transparentTarget, { label: 'helpers' });
      if (helper.debug.isEnabled) {
        if (b.drawing) helper.debug.syncVisibility();
        for (const scene of helper.debug.scenes) {
          this.drawPass.renderBlended(b, renderer, scene, camera);
        }
      }
      if (helper.pointer.isEnabled) {
        if (b.drawing) {
          helper.pointer.setCamera(camera);
          renderer.update(helper.pointer.camera, helper.pointer.scene, frame);
        }
        this.drawPass.renderBlended(b, renderer, helper.pointer.scene, helper.pointer.camera);
      }
      this.drawPass.depthTextureOpaque.detachFramebuffer(this.transparentTarget.framebuffer, 'depth');
      b.begin(this.transparentTarget, { label: 'helpers' });
    }
    if (helper.handle.isEnabled) {
      this.drawPass.renderBlended(b, renderer, helper.handle.scene, camera);
    }
    if (helper.camera.isEnabled) {
      if (b.drawing) {
        helper.camera.update(camera);
        renderer.update(helper.camera.camera, helper.camera.scene, frame);
      }
      this.drawPass.renderBlended(b, renderer, helper.camera.scene, helper.camera.camera);
    }
  }

  shouldRender(props: IlluminationProps) {
    return this._supported && props.enabled && this._iteration < this.getMaxIterations(props);
  }

  setSize(width: number, height: number) {
    if (!this._supported) return;

    const w = this.outputTarget.getWidth();
    const h = this.outputTarget.getHeight();

    if (width !== w || height !== h) {
      this.tracing.setSize(width, height);

      this.transparentTarget.setSize(width, height);
      this.outputTarget.setSize(width, height);
      this.denoiser?.output.define(width, height);

      ValueCell.update(
        this.copyRenderable.values.uTexSize,
        Vec2.set(this.copyRenderable.values.uTexSize.ref.value, width, height),
      );
      ValueCell.update(
        this.composeRenderable.values.uTexSize,
        Vec2.set(this.composeRenderable.values.uTexSize.ref.value, width, height),
      );

      this.multiSampleComposeTarget.setSize(width, height);
      this.multiSampleHoldTarget.setSize(width, height);
      ValueCell.update(
        this.multiSampleCompose.values.uTexSize,
        Vec2.set(this.multiSampleCompose.values.uTexSize.ref.value, width, height),
      );
    }

    this.drawPass.setSize(width, height);
  }

  restart(clearAdjustedProps = false) {
    if (!this._supported) return;

    this.tracing.restart(clearAdjustedProps);
    this._iteration = 0;
    this.prevSampleIndex = -1;
  }

  private renderInternal(
    b: PassBuilder,
    ctx: RenderContext<Camera>,
    props: Props,
    toDrawingBuffer: boolean,
    forceRenderInput: boolean,
    iteration: number,
  ) {
    if (b.drawing && !this.shouldRender(props.illumination)) return;

    const timing = isTimingMode && b.drawing;
    if (timing) {
      this.graphics.timer.group('IlluminationPass.render', {
        note: `iteration ${iteration + 1} of ${this.getMaxIterations(props.illumination)}`,
      });
    }

    this.tracing.render(b, ctx, props.transparentBackground, props.illumination, iteration, forceRenderInput);
    this.compose(b, ctx, props, toDrawingBuffer, iteration === 0 || forceRenderInput);

    if (!b.drawing) return;
    this._iteration += 1;
    if (timing) this.graphics.timer.groupEnd('IlluminationPass.render');

    this.graphics.flush();
  }

  private setupRenderer(ctx: RenderContext<Camera>, props: Props) {
    const { renderer, camera, scene, frame } = ctx;
    renderer.setTransparentBackground(props.transparentBackground);
    renderer.setDrawingBufferSize(this.tracing.composeTarget.getWidth(), this.tracing.composeTarget.getHeight());
    renderer.setPixelRatio(this.graphics.pixelRatio);
    const { x, y, width, height } = camera.viewport;
    renderer.setViewport(x, y, width, height);
    renderer.update(camera, scene, frame);
  }

  private getDenoiser() {
    if (!this.denoiser) {
      const { graphics, tracing } = this;
      const { accumulateTarget } = tracing;
      const output = this.drawPass.targets.texture('image-float16', 'rgba', 'fp16', 'nearest');
      output.define(accumulateTarget.getWidth(), accumulateTarget.getHeight());
      this.denoiser = new Denoiser(graphics.compute!, accumulateTarget.texture, tracing.normalTextureOpaque, output);
    }
    return this.denoiser;
  }

  private compose(
    b: PassBuilder,
    ctx: RenderContext<Camera>,
    props: Props,
    toDrawingBuffer: boolean,
    forceRenderInput: boolean,
  ) {
    const { renderer, camera, scene, frame } = ctx;
    const { x, y, width, height } = camera.viewport;

    if (forceRenderInput) {
      // render color & depth
      if (b.drawing) this.setupRenderer(ctx, props);
      else renderer.setTransparentBackground(props.transparentBackground);
      this.renderInput(b, renderer, camera, scene, props);
      this.renderHelpers(b, ctx);
    }

    if (b.drawing) {
      this.graphics.setRenderState(RenderState.Default);
      this.graphics.setViewport(x, y, width, height);
      this.graphics.setScissor(x, y, width, height);
    }

    const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

    const antialiasingEnabled = AntialiasingPass.isEnabled(props.postprocessing);
    const outlinesEnabled = OutlinePass.isEnabled(props.postprocessing) && !props.illumination.ignoreOutline;
    const occlusionEnabled = PostprocessingPass.isTransparentSsaoEnabled(scene, props.postprocessing);
    const bloomEnabled = BloomPass.isActive(scene, props.postprocessing);
    const dofEnabled = DofPass.isEnabled(props.postprocessing);

    const hasTransparent = scene.opacityAverage < 1 || scene.volumes.renderables.length > 0;

    let needsUpdateCompose = false;

    if (this.composeRenderable.values.dOutlineEnable.ref.value !== outlinesEnabled) {
      needsUpdateCompose = true;
      ValueCell.update(this.composeRenderable.values.dOutlineEnable, outlinesEnabled);
    }

    if (outlinesEnabled && props.postprocessing.outline.name === 'on') {
      const { transparentOutline, outlineScale } = this.drawPass.postprocessing.outline.update(
        camera,
        props.postprocessing.outline.params,
        this.drawPass.depthTargetTransparent.texture,
        this.drawPass.depthTextureOpaque,
      );
      this.drawPass.postprocessing.outline.render(b);

      ValueCell.update(
        this.composeRenderable.values.uOutlineColor,
        Color.toVec3Normalized(
          this.composeRenderable.values.uOutlineColor.ref.value,
          props.postprocessing.outline.params.color,
        ),
      );

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

    if (occlusionEnabled) {
      const [width, height] = this.composeRenderable.values.uTexSize.ref.value;
      const ssaoTexture = this.drawPass.postprocessing.ssao.ssaoDepthTexture;
      if (updateOcclusionUpsample(this.composeRenderable.values, ssaoTexture, width, height)) needsUpdateCompose = true;
    }

    if (occlusionEnabled && props.postprocessing.occlusion.name === 'on') {
      ValueCell.update(
        this.composeRenderable.values.uOcclusionColor,
        Color.toVec3Normalized(
          this.composeRenderable.values.uOcclusionColor.ref.value,
          props.postprocessing.occlusion.params.color,
        ),
      );
    }

    // helpers are drawn into the transparent target, so it must be blended even without transparent geometry
    const blendTransparency = hasTransparent || this.hasHelpers(ctx.helper);
    if (this.composeRenderable.values.dBlendTransparency.ref.value !== blendTransparency) {
      needsUpdateCompose = true;
      ValueCell.update(this.composeRenderable.values.dBlendTransparency, blendTransparency);
    }

    ValueCell.updateIfChanged(this.composeRenderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.composeRenderable.values.uFar, camera.far);
    ValueCell.updateIfChanged(this.composeRenderable.values.uFogFar, camera.fogFar);
    ValueCell.updateIfChanged(this.composeRenderable.values.uFogNear, camera.fogNear);
    ValueCell.update(
      this.composeRenderable.values.uFogColor,
      Color.toVec3Normalized(this.composeRenderable.values.uFogColor.ref.value, renderer.props.backgroundColor),
    );
    if (this.composeRenderable.values.dOrthographic.ref.value !== orthographic) {
      ValueCell.update(this.composeRenderable.values.dOrthographic, orthographic);
      needsUpdateCompose = true;
    }

    // bloom composite — compose applies it inline below

    if (b.drawing) {
      renderer.setDrawingBufferSize(this.tracing.composeTarget.getWidth(), this.tracing.composeTarget.getHeight());
      renderer.setPixelRatio(this.graphics.pixelRatio);
      renderer.setViewport(x, y, width, height);
      renderer.update(camera, scene, frame);
    }

    let bloomActive = false;
    if (bloomEnabled && props.postprocessing.bloom.name === 'on') {
      const params = props.postprocessing.bloom.params;
      const bloom = this.drawPass.postprocessing.bloom;
      const emissiveBloom = params.mode === 'emissive';
      const skip = emissiveBloom && scene.emissiveAverage === 0;

      if (!skip) {
        if (emissiveBloom) {
          this.drawPass.renderEmissiveBloom(b, renderer, camera, scene, params.transparency);
        }
        // transparent color is already cleared in renderInput and may hold helpers, so only reset depth
        if (scene.opacityAverage >= 1 && b.drawing) {
          b.begin(this.drawPass.depthTargetTransparent, { label: 'transparent depth clear' });
          renderer.clearDepth(true);
        }
        // seed luminosity bloom from the shaded (lit) color, like the standard path's colorTarget
        bloom.update(
          this.tracing.shadedTextureOpaque,
          this.transparentTarget.texture,
          bloom.emissiveTarget.texture,
          this.drawPass.depthTextureOpaque,
          this.drawPass.depthTextureTransparent,
          params,
          camera,
          false,
          // the transparent target may hold helpers, see `renderInput`
          true,
        );
        bloom.render(b, camera.viewport);
        bloomActive = true;
      }
    }
    if (this.composeRenderable.values.dBloomEnable.ref.value !== bloomActive) {
      needsUpdateCompose = true;
      ValueCell.update(this.composeRenderable.values.dBloomEnable, bloomActive);
    }

    // denoise

    const denoiseThreshold =
      props.multiSample.mode === 'on'
        ? props.illumination.denoiseThreshold[0]
        : lerp(
            props.illumination.denoiseThreshold[1],
            props.illumination.denoiseThreshold[0],
            clamp(this.iteration / (this.getMaxIterations(props.illumination) / 2), 0, 1),
          );
    // a kernel where compute shaders run, a texel read and its normal decoded once per tile instead of by each of its
    // 169 neighbours; the compose program reads its output
    const computeDenoise = props.illumination.denoise && isComputeDenoiseSupported(this.graphics);
    if (computeDenoise && b.drawing) this.getDenoiser().run(denoiseThreshold);
    const composeColor = computeDenoise && this.denoiser ? this.denoiser.output : this.tracing.accumulateTarget.texture;
    if (this.composeRenderable.values.tColor.ref.value !== composeColor) {
      ValueCell.update(this.composeRenderable.values.tColor, composeColor);
      needsUpdateCompose = true;
    }

    // background

    // the image is always kept offscreen, so that marking can be redrawn over it without restarting
    b.begin(this.tracing.composeTarget, { label: 'compose' });
    this._colorTarget = this.tracing.composeTarget;

    const { background } = this.drawPass.postprocessing;
    if (b.drawing) {
      background.update(camera, props.postprocessing.background);
      background.clear(props.postprocessing.background, props.transparentBackground, renderer.props.backgroundColor);
      background.render(b, props.postprocessing.background);
    } else {
      // the props may be hypothetical (`Canvas3D.prepare`): the live background stays
      background.prepare(b, props.postprocessing.background);
    }

    // compose

    // premultiplied scene composited over the background (cleared color or environment) with
    // premultiplied-over — keeps transparency weight (vs SRC_ALPHA's color*a^2) and lets bloom occlude
    const backgroundVariant = this.drawPass.postprocessing.background.isEnabled(props.postprocessing);
    const compositeBackground = backgroundVariant || (bloomActive && !props.transparentBackground);
    ValueCell.updateIfChanged(
      this.composeRenderable.values.uTransparentBackground,
      props.transparentBackground || compositeBackground,
    );
    const programDenoise = props.illumination.denoise && !computeDenoise;
    if (this.composeRenderable.values.dDenoise.ref.value !== programDenoise) {
      ValueCell.update(this.composeRenderable.values.dDenoise, programDenoise);
      needsUpdateCompose = true;
    }
    ValueCell.updateIfChanged(this.composeRenderable.values.uDenoiseThreshold, denoiseThreshold);
    if (needsUpdateCompose) this.composeRenderable.update();
    b.quad(this.composeRenderable, compositeBackground ? RenderState.Premultiplied : RenderState.Default);

    //

    let swapTarget = this.outputTarget;

    if (antialiasingEnabled) {
      this.drawPass.antialiasing.render(
        b,
        camera,
        this.tracing.composeTarget.texture,
        this.outputTarget,
        props.postprocessing,
      );
      this._colorTarget = this.outputTarget;
      swapTarget = this.tracing.composeTarget;
    }

    if (dofEnabled && props.postprocessing.dof.name === 'on') {
      this.drawPass.dof.update(
        camera,
        this._colorTarget.texture,
        this.drawPass.depthTextureOpaque,
        this.drawPass.depthTargetTransparent.texture,
        props.postprocessing.dof.params,
        scene.boundingSphereVisible,
      );
      this.drawPass.dof.render(b, camera.viewport, swapTarget);
      this._colorTarget = swapTarget;
    }

    //

    // in 'on' mode, marking is blended over the accumulated samples instead
    if (props.multiSample.mode === 'on' && !toDrawingBuffer) return;

    this.drawPass.marking.present(b, ctx, props, {
      base: this._colorTarget,
      toDrawingBuffer,
      offsets: SingleSample,
      restart: forceRenderInput,
      samples: 1,
      shading: this.markingShading,
    });
  }

  private prevSampleIndex = -1;

  private renderMultiSample(
    b: PassBuilder,
    ctx: RenderContext<Camera>,
    props: Props,
    toDrawingBuffer: boolean,
    currentIteration: number,
  ) {
    const { camera } = ctx;
    const { multiSampleCompose, multiSampleComposeTarget, multiSampleHoldTarget, graphics } = this;

    // based on the Multisample Anti-Aliasing Render Pass
    // contributed to three.js by bhouston / http://clara.io/
    //
    // This manual approach to MSAA re-renders the scene once for
    // each sample with camera jitter and accumulates the results.
    const offsetList = getJitterOffsets(props.multiSample.sampleLevel);

    const maxIterations = this.getMaxIterations(props.illumination);
    const iteration = Math.min(currentIteration, maxIterations);

    const sampleIndex = Math.floor((iteration / maxIterations) * offsetList.length);

    const timing = isTimingMode && b.drawing;
    if (timing) {
      graphics.timer.group('IlluminationPass.renderMultiSample', {
        note: `sampleIndex ${sampleIndex + 1} of ${offsetList.length}`,
      });
    }

    const { x, y, width, height } = camera.viewport;
    const sampleWeight = 1.0 / maxIterations;

    // the compose target holds `w * (s_1 + ... + s_k) + (1 - k * w) * s_0` after iteration `k` (samples `s`, weight
    // `w`): the first sample, the hold, until every following one has its weight; kept by adding `w * (s_k - s_0)`
    // per iteration (the first sample has its weight in the sum of the last iteration, `k * w = 1`)
    if (iteration === 0) {
      this.renderInternal(b, ctx, props, false, true, iteration);
      ValueCell.update(multiSampleCompose.values.uWeight, 1.0);
      ValueCell.update(multiSampleCompose.values.tColor, this._colorTarget.texture);
      multiSampleCompose.update();

      b.begin(multiSampleHoldTarget, { label: 'multi-sample hold' });
      if (b.drawing) {
        graphics.setRenderState(RenderState.Default);
        graphics.setViewport(x, y, width, height);
        graphics.setScissor(x, y, width, height);
      }
      b.quad(multiSampleCompose, RenderState.Default);

      b.begin(multiSampleComposeTarget, { label: 'multi-sample compose' });
      if (b.drawing) {
        graphics.setViewport(x, y, width, height);
        graphics.setScissor(x, y, width, height);
      }
      b.quad(multiSampleCompose, RenderState.Default);
    } else {
      // render the scene multiple times, each slightly jitter offset
      // from the last and accumulate the results.
      if (b.drawing) setJitter(camera, offsetList[sampleIndex]);

      // render scene
      this.renderInternal(b, ctx, props, false, !b.drawing || this.prevSampleIndex !== sampleIndex, iteration);

      // add the sample's weight and take it from the first sample (a float target, not clamped)
      b.begin(multiSampleComposeTarget, { label: 'multi-sample compose' });
      if (b.drawing) {
        graphics.setViewport(x, y, width, height);
        graphics.setScissor(x, y, width, height);
      }
      ValueCell.update(multiSampleCompose.values.tColor, this._colorTarget.texture);
      ValueCell.update(multiSampleCompose.values.uWeight, sampleWeight);
      multiSampleCompose.update();
      b.quad(multiSampleCompose, RenderState.Additive);
      ValueCell.update(multiSampleCompose.values.tColor, multiSampleHoldTarget.texture);
      ValueCell.update(multiSampleCompose.values.uWeight, -sampleWeight);
      multiSampleCompose.update();
      b.quad(multiSampleCompose, RenderState.Additive);
    }

    if (b.drawing) this.prevSampleIndex = sampleIndex;

    if (!toDrawingBuffer) {
      b.begin(this.colorTarget, { label: 'multi-sample output' });
      if (this.copyRenderable.values.tColor.ref.value !== multiSampleComposeTarget.texture) {
        ValueCell.update(this.copyRenderable.values.tColor, multiSampleComposeTarget.texture);
        this.copyRenderable.update();
      }
      b.quad(this.copyRenderable, RenderState.Default);
    }

    if (b.drawing) clearJitter(camera);

    const base = toDrawingBuffer ? multiSampleComposeTarget : this._colorTarget;
    // a single sample on the first iteration, which is also the one rendered while the camera moves
    const samples = iteration === 0 ? 1 : Math.ceil(offsetList.length / maxIterations);
    this.drawPass.marking.present(b, ctx, props, {
      base,
      toDrawingBuffer,
      offsets: offsetList,
      restart: iteration === 0,
      samples,
      shading: this.markingShading,
    });
    if (timing) graphics.timer.groupEnd('IlluminationPass.renderMultiSample');
  }

  render(ctx: RenderContext<Camera>, props: Props, toDrawingBuffer: boolean) {
    this.plan(this.renderBuilder, ctx, props, toDrawingBuffer);
  }

  /**
   * What `render` draws, see `PassBuilder`; preparing covers the first and the following iterations (the input is
   * drawn on the first and when a multi-sample offset changes) and keeps the iteration and the output target
   */
  plan(b: PassBuilder, ctx: RenderContext<Camera>, props: Props, toDrawingBuffer: boolean) {
    if (!this._supported) return;

    if (b.drawing) {
      if (props.multiSample.mode === 'on') {
        this.renderMultiSample(b, ctx, props, toDrawingBuffer, this._iteration);
      } else {
        this.renderInternal(b, ctx, props, toDrawingBuffer, false, this._iteration);
      }
      return;
    }

    const colorTarget = this._colorTarget;
    for (const iteration of [0, 1]) {
      if (props.multiSample.mode === 'on') {
        this.renderMultiSample(b, ctx, props, toDrawingBuffer, iteration);
      } else {
        this.renderInternal(b, ctx, props, toDrawingBuffer, iteration === 0, iteration);
      }
    }
    this._colorTarget = colorTarget;
  }
}

//

export const ComposeSchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tNormal: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tShaded: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tTransparentColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dBlendTransparency: DefineSpec('boolean'),
  tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tSsaoDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  /** the SSAO targets have a lower resolution, see `getSsaoUpsampled` */
  dOcclusionUpsample: DefineSpec('boolean'),
  uOcclusionTexSize: UniformSpec('v2'),
  tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tOutlines: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tBloom: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  dBloomEnable: DefineSpec('boolean'),
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
export const ComposeShaderCode = ShaderCode('illumination-compose');
type ComposeRenderable = ComputeRenderable<Values<typeof ComposeSchema>>;

function getComposeRenderable(
  ctx: GraphicsContext,
  colorTexture: Texture,
  normalTexture: Texture,
  shadedTexture: Texture,
  depthTextureOpaque: Texture,
  depthTextureTransparent: Texture,
  outlinesTexture: Texture,
  transparentColorTexture: Texture,
  ssaoDepthOpaqueTexture: Texture,
  ssaoDepthTransparentTexture: Texture,
  bloomTexture: Texture,
  transparentOutline: boolean,
): ComposeRenderable {
  const values: Values<typeof ComposeSchema> = {
    ...QuadValues,
    tColor: ValueCell.create(colorTexture),
    tNormal: ValueCell.create(normalTexture),
    tShaded: ValueCell.create(shadedTexture),
    tTransparentColor: ValueCell.create(transparentColorTexture),
    dBlendTransparency: ValueCell.create(true),
    tSsaoDepth: ValueCell.create(ssaoDepthOpaqueTexture),
    tSsaoDepthTransparent: ValueCell.create(ssaoDepthTransparentTexture),
    dOcclusionUpsample: ValueCell.create(false),
    uOcclusionTexSize: ValueCell.create(Vec2()),
    tDepthOpaque: ValueCell.create(depthTextureOpaque),
    tDepthTransparent: ValueCell.create(depthTextureTransparent),
    tOutlines: ValueCell.create(outlinesTexture),
    tBloom: ValueCell.create(bloomTexture),
    dBloomEnable: ValueCell.create(false),
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
  const renderItem = ctx.createComputeRenderItem('triangles', ComposeShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

export const MultiSampleComposeSchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),
  uWeight: UniformSpec('f'),
};
export const MultiSampleComposeShaderCode = ShaderCode('illumination-multi-sample-compose');
type MultiSampleComposeRenderable = ComputeRenderable<Values<typeof MultiSampleComposeSchema>>;

function getMultiSampleComposeRenderable(ctx: GraphicsContext, colorTexture: Texture): MultiSampleComposeRenderable {
  const values: Values<typeof MultiSampleComposeSchema> = {
    ...QuadValues,
    tColor: ValueCell.create(colorTexture),
    uTexSize: ValueCell.create(Vec2.create(colorTexture.getWidth(), colorTexture.getHeight())),
    uWeight: ValueCell.create(1.0),
  };

  const schema = { ...MultiSampleComposeSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', MultiSampleComposeShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
