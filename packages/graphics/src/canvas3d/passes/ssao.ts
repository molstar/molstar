/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PassBuilder } from './pass-builder.js';
import { type CopyRenderable, createCopyRenderable, QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { TextureSpec, type Values, UniformSpec, DefineSpec } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import { deepEqual, ValueCell } from '@molstar/core/util';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { ICamera } from '../camera.js';
import type { Scene } from '@molstar/graphics/gl/scene';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import { Color } from '@molstar/core/util/color';
import { isTimingMode } from '@molstar/core/util/debug';
import type { PostprocessingProps } from './postprocessing.js';
import { PCG } from '@molstar/core/data/util/hash-functions';
import type { PassTargets } from './pass-targets.js';

/** max. array lengths in uniform blocks, see `SsaoParams` */
const MaxSamples = 256;
const MaxLevels = 16;
const MaxBlurKernelSize = 25;

export const SsaoParams = {
  samples: PD.Numeric(32, { min: 1, max: MaxSamples, step: 1 }),
  multiScale: PD.MappedStatic(
    'off',
    {
      on: PD.Group({
        levels: PD.ObjectList(
          {
            radius: PD.Numeric(5, { min: 0, max: 20, step: 0.1 }, { description: 'Final occlusion radius is 2^x' }),
            bias: PD.Numeric(1, { min: 0, max: 3, step: 0.1 }),
          },
          (o) => `${o.radius}, ${o.bias}`,
          {
            defaultValue: [
              { radius: 2, bias: 1 },
              { radius: 5, bias: 1 },
              { radius: 8, bias: 1 },
              { radius: 11, bias: 1 },
            ],
          },
        ),
        nearThreshold: PD.Numeric(10, { min: 0, max: 50, step: 1 }),
        farThreshold: PD.Numeric(1500, { min: 0, max: 10000, step: 100 }),
      }),
      off: PD.Group({}),
    },
    { cycle: true },
  ),
  radius: PD.Numeric(
    5,
    { min: 0, max: 20, step: 0.1 },
    { description: 'Final occlusion radius is 2^x', hideIf: (p) => p?.multiScale.name === 'on' },
  ),
  bias: PD.Numeric(0.8, { min: 0, max: 3, step: 0.1 }),
  blurKernelSize: PD.Numeric(15, { min: 1, max: MaxBlurKernelSize, step: 2 }),
  blurDepthBias: PD.Numeric(0.5, { min: 0, max: 1, step: 0.01 }),
  resolutionScale: PD.Numeric(
    0.5,
    { min: 0.1, max: 1, step: 0.05 },
    {
      description:
        'Resolution of the occlusion relative to the drawing buffer, at most that of CSS pixels (high pixel ratios), ' +
        'upsampled with the depth',
    },
  ),
  color: PD.Color(Color(0x000000)),
  transparentThreshold: PD.Numeric(0.4, { min: 0, max: 1, step: 0.05 }),
};

export type SsaoProps = PD.Values<typeof SsaoParams>;
const DefaultResolutionScale = SsaoParams.resolutionScale.defaultValue;

/** Values of passes reading the SSAO targets, which are upsampled depth-aware below the full resolution */
export type OcclusionUpsampleValues = {
  dOcclusionUpsample: ValueCell<boolean>;
  uOcclusionTexSize: ValueCell<Vec2>;
};

/** Updates `values` for SSAO targets of `ssaoTexture`'s size read at `width` x `height`; true if the define changed */
export function updateOcclusionUpsample(
  values: OcclusionUpsampleValues,
  ssaoTexture: Texture,
  width: number,
  height: number,
) {
  const w = ssaoTexture.getWidth();
  const h = ssaoTexture.getHeight();
  const size = values.uOcclusionTexSize.ref.value;
  if (size[0] !== w || size[1] !== h) ValueCell.update(values.uOcclusionTexSize, Vec2.set(size, w, h));
  const upsample = w !== width || h !== height;
  if (values.dOcclusionUpsample.ref.value === upsample) return false;
  ValueCell.update(values.dOcclusionUpsample, upsample);
  return true;
}

type Levels = {
  count: number;
  radius: number[];
  bias: number[];
};

function getLevels(props: { radius: number; bias: number }[], scale: number, levels?: Levels): Levels {
  const count = Math.min(props.length, MaxLevels);
  const { radius, bias } = levels || {
    radius: new Array(count * 3).fill(0),
    bias: new Array(count * 3).fill(0),
  };
  props = props.slice().sort((a, b) => a.radius - b.radius);
  for (let i = 0; i < count; ++i) {
    const p = props[i];
    radius[i] = Math.pow(2, p.radius) * scale;
    bias[i] = p.bias;
  }
  return { count, radius, bias };
}

export class SsaoPass {
  private readonly graphics: GraphicsContext;
  static isEnabled(props: PostprocessingProps) {
    return props.enabled && props.occlusion.name !== 'off';
  }

  static isTransparentEnabled(scene: Scene, props: SsaoProps) {
    return scene.opacityAverage < 1 && scene.transparencyMin < props.transparentThreshold;
  }

  private readonly framebuffer: Framebuffer;
  private readonly blurFirstPassFramebuffer: Framebuffer;
  private readonly blurSecondPassFramebuffer: Framebuffer;

  /** no depth textures, the opaque depth itself is packed */
  private readonly packedDepth: boolean;
  /** half and quarter depth targets are packed into RGBA8, without float render targets */
  private readonly packedDepthLevels: boolean;

  private readonly depthHalfTargetOpaque: RenderTarget;
  private readonly depthHalfRenderableOpaque: CopyRenderable;

  private readonly depthQuarterTargetOpaque: RenderTarget;
  private readonly depthQuarterRenderableOpaque: CopyRenderable;

  private readonly depthHalfTargetTransparent: RenderTarget;
  private readonly depthHalfRenderableTransparent: CopyRenderable;

  private readonly depthQuarterTargetTransparent: RenderTarget;
  private readonly depthQuarterRenderableTransparent: CopyRenderable;

  readonly ssaoDepthTexture: Texture;
  readonly ssaoDepthTransparentTexture: Texture;

  private readonly depthBlurProxyTexture: Texture;

  private depthTextureOpaque: Texture;
  private depthTextureTransparent: Texture;

  private readonly renderable: SsaoRenderable;
  private readonly blurFirstPassRenderable: SsaoBlurRenderable;
  private readonly blurSecondPassRenderable: SsaoBlurRenderable;

  private nSamples: number;
  private blurKernelSize: number;
  /** a depth bias of 0 blurs nothing (no sample is within it): the blur passes are skipped */
  private blur = true;
  private texSize: [number, number];
  private invProjection = Mat4.identity();

  private ssaoScale: number;
  /** `SsaoProps.resolutionScale` of the last update */
  private resolutionScale = DefaultResolutionScale;
  private calcSsaoScale(resolutionScale: number) {
    // at most the resolution of CSS pixels for high pixel ratios
    return Math.min(1, resolutionScale, 1 / this.graphics.pixelRatio);
  }

  private levelsCameraScale = -1;
  private levels: { radius: number; bias: number }[];

  // the SSAO pass reads the full-resolution depth at any SSAO resolution: a downsampled copy (nearest texels) widens
  // silhouettes, more samples near them would hit the object and the occlusion would grow with lower resolutions
  private getDepthTexture() {
    return this.depthTextureOpaque;
  }

  private getTransparentDepthTexture() {
    return this.depthTextureTransparent;
  }

  private updatePackedDepth() {
    ValueCell.updateIfChanged(this.renderable.values.dPackedDepth, this.packedDepth);
    ValueCell.updateIfChanged(
      this.depthHalfRenderableOpaque.values.dPackDepth,
      this.packedDepthLevels && !this.packedDepth,
    );
  }

  constructor(
    targets: PassTargets,
    width: number,
    height: number,
    packedDepth: boolean,
    depthTextureOpaque: Texture,
    depthTextureTransparent: Texture,
  ) {
    const graphics = (this.graphics = targets.graphics);
    const { floatRenderTarget, floatTextureLinear } = graphics.capabilities;

    this.packedDepth = packedDepth;
    this.packedDepthLevels = packedDepth || !floatRenderTarget;

    this.depthTextureOpaque = depthTextureOpaque;
    this.depthTextureTransparent = depthTextureTransparent;

    this.nSamples = 1;
    this.blurKernelSize = 1;
    this.ssaoScale = this.calcSsaoScale(this.resolutionScale);
    this.texSize = [width, height];
    this.levels = [];

    this.framebuffer = graphics.resources.framebuffer();
    this.blurFirstPassFramebuffer = graphics.resources.framebuffer();
    this.blurSecondPassFramebuffer = graphics.resources.framebuffer();

    const sw = Math.floor(width * this.ssaoScale);
    const sh = Math.floor(height * this.ssaoScale);

    const hw = Math.max(1, Math.floor(sw * 0.5));
    const hh = Math.max(1, Math.floor(sh * 0.5));

    const qw = Math.max(1, Math.floor(sw * 0.25));
    const qh = Math.max(1, Math.floor(sh * 0.25));

    const filter = floatTextureLinear ? 'linear' : 'nearest';
    const createDepthTarget = (w: number, h: number) =>
      this.packedDepthLevels
        ? targets.createRenderTarget(w, h, 'none', 'uint8', 'nearest', 'rgba')
        : targets.createRenderTarget(w, h, 'none', 'float32', filter, 'r');
    const depthTexture = this.getDepthTexture();
    this.depthHalfTargetOpaque = createDepthTarget(hw, hh);
    this.depthHalfRenderableOpaque = createCopyRenderable(graphics, depthTexture);

    this.depthQuarterTargetOpaque = createDepthTarget(qw, qh);
    this.depthQuarterRenderableOpaque = createCopyRenderable(graphics, this.depthHalfTargetOpaque.texture);

    const transparentDepthTexture = this.getTransparentDepthTexture();
    this.depthHalfTargetTransparent = targets.createRenderTarget(hw, hh, 'none', 'uint8', 'linear', 'rgba');
    this.depthHalfRenderableTransparent = createCopyRenderable(graphics, transparentDepthTexture);

    this.depthQuarterTargetTransparent = targets.createRenderTarget(qw, qh, 'none', 'uint8', 'linear', 'rgba');
    this.depthQuarterRenderableTransparent = createCopyRenderable(graphics, this.depthHalfTargetTransparent.texture);

    this.ssaoDepthTexture = targets.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    this.ssaoDepthTexture.define(sw, sh);
    this.ssaoDepthTexture.attachFramebuffer(this.framebuffer, 'color0');

    this.ssaoDepthTransparentTexture = targets.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    this.ssaoDepthTransparentTexture.define(sw, sh);

    this.depthBlurProxyTexture = targets.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    this.depthBlurProxyTexture.define(sw, sh);
    this.depthBlurProxyTexture.attachFramebuffer(this.blurFirstPassFramebuffer, 'color0');

    this.renderable = getSsaoRenderable(
      graphics,
      depthTexture,
      this.depthHalfTargetOpaque.texture,
      this.depthQuarterTargetOpaque.texture,
      transparentDepthTexture,
      this.depthHalfTargetTransparent.texture,
      this.depthQuarterTargetTransparent.texture,
      this.packedDepthLevels,
    );
    this.blurFirstPassRenderable = getSsaoBlurRenderable(graphics, this.ssaoDepthTransparentTexture, 'horizontal');
    this.blurSecondPassRenderable = getSsaoBlurRenderable(graphics, this.depthBlurProxyTexture, 'vertical');

    this.updatePackedDepth();
    this.renderable.update();
    this.depthHalfRenderableOpaque.update();

    // the sizes of the level copies, the passes and the blur as after any resize (the renderables above take the sizes
    // of their inputs)
    this.texSize = [0, 0];
    this.setSize(width, height);
  }

  getByteCount() {
    return (
      this.depthHalfTargetOpaque.getByteCount() +
      this.depthQuarterTargetOpaque.getByteCount() +
      this.depthHalfTargetTransparent.getByteCount() +
      this.depthQuarterTargetTransparent.getByteCount() +
      this.ssaoDepthTexture.getByteCount() +
      this.ssaoDepthTransparentTexture.getByteCount() +
      this.depthBlurProxyTexture.getByteCount()
    );
  }

  setSize(width: number, height: number) {
    const [w, h] = this.texSize;
    const ssaoScale = this.calcSsaoScale(this.resolutionScale);
    if (width !== w || height !== h || this.ssaoScale !== ssaoScale) {
      this.texSize.splice(0, 2, width, height);
      this.ssaoScale = ssaoScale;

      const sw = Math.floor(width * this.ssaoScale);
      const sh = Math.floor(height * this.ssaoScale);
      this.ssaoDepthTexture.define(sw, sh);
      this.ssaoDepthTransparentTexture.define(sw, sh);
      this.depthBlurProxyTexture.define(sw, sh);

      const hw = Math.max(1, Math.floor(sw * 0.5));
      const hh = Math.max(1, Math.floor(sh * 0.5));
      const qw = Math.max(1, Math.floor(sw * 0.25));
      const qh = Math.max(1, Math.floor(sh * 0.25));

      this.depthHalfTargetOpaque.setSize(hw, hh);
      this.depthQuarterTargetOpaque.setSize(qw, qh);

      ValueCell.update(
        this.depthHalfRenderableOpaque.values.uTexSize,
        Vec2.set(this.depthHalfRenderableOpaque.values.uTexSize.ref.value, hw, hh),
      );
      ValueCell.update(
        this.depthQuarterRenderableOpaque.values.uTexSize,
        Vec2.set(this.depthQuarterRenderableOpaque.values.uTexSize.ref.value, qw, qh),
      );

      this.depthHalfTargetTransparent.setSize(hw, hh);
      this.depthQuarterTargetTransparent.setSize(qw, qh);

      ValueCell.update(
        this.depthHalfRenderableTransparent.values.uTexSize,
        Vec2.set(this.depthHalfRenderableTransparent.values.uTexSize.ref.value, hw, hh),
      );
      ValueCell.update(
        this.depthQuarterRenderableTransparent.values.uTexSize,
        Vec2.set(this.depthQuarterRenderableTransparent.values.uTexSize.ref.value, qw, qh),
      );

      ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, sw, sh));
      ValueCell.update(
        this.blurFirstPassRenderable.values.uTexSize,
        Vec2.set(this.blurFirstPassRenderable.values.uTexSize.ref.value, sw, sh),
      );
      ValueCell.update(
        this.blurSecondPassRenderable.values.uTexSize,
        Vec2.set(this.blurSecondPassRenderable.values.uTexSize.ref.value, sw, sh),
      );

      const depthTexture = this.getDepthTexture();
      const transparentDepthTexture = this.getTransparentDepthTexture();

      ValueCell.update(this.depthHalfRenderableOpaque.values.tColor, depthTexture);
      ValueCell.update(this.depthHalfRenderableTransparent.values.tColor, transparentDepthTexture);

      ValueCell.update(this.renderable.values.tDepth, depthTexture);
      ValueCell.update(this.renderable.values.tDepthTransparent, transparentDepthTexture);

      this.depthHalfRenderableOpaque.update();
      this.depthHalfRenderableTransparent.update();
      this.renderable.update();
    }
  }

  update(camera: ICamera, scene: Scene, props: SsaoProps, illuminationMode = false) {
    let needsUpdateSsao = false;
    let needsUpdateSsaoBlur = false;
    let needsUpdateDepthHalf = false;

    const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;
    const invProjection = Mat4.invert(this.invProjection, camera.projection);

    const [w, h] = this.texSize;
    const v = camera.viewport;

    ValueCell.update(this.renderable.values.uProjection, camera.projection);
    ValueCell.update(this.renderable.values.uInvProjection, invProjection);
    ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);

    const b = this.renderable.values.uBounds;
    const s = this.ssaoScale;
    Vec4.set(
      b.ref.value,
      Math.floor(v.x * s) / (w * s),
      Math.floor(v.y * s) / (h * s),
      Math.ceil((v.x + v.width) * s) / (w * s),
      Math.ceil((v.y + v.height) * s) / (h * s),
    );
    ValueCell.update(b, b.ref.value);
    ValueCell.update(this.blurFirstPassRenderable.values.uBounds, b.ref.value);
    ValueCell.update(this.blurSecondPassRenderable.values.uBounds, b.ref.value);

    ValueCell.updateIfChanged(this.blurFirstPassRenderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.blurSecondPassRenderable.values.uNear, camera.near);

    ValueCell.updateIfChanged(this.blurFirstPassRenderable.values.uFar, camera.far);
    ValueCell.updateIfChanged(this.blurSecondPassRenderable.values.uFar, camera.far);

    ValueCell.update(this.blurFirstPassRenderable.values.uInvProjection, invProjection);
    ValueCell.update(this.blurSecondPassRenderable.values.uInvProjection, invProjection);

    ValueCell.update(this.blurFirstPassRenderable.values.uBlurDepthBias, props.blurDepthBias * camera.scale);
    ValueCell.update(this.blurSecondPassRenderable.values.uBlurDepthBias, props.blurDepthBias * camera.scale);
    this.blur = props.blurDepthBias > 0;

    if (this.blurFirstPassRenderable.values.dOrthographic.ref.value !== orthographic) {
      needsUpdateSsaoBlur = true;
      ValueCell.update(this.blurFirstPassRenderable.values.dOrthographic, orthographic);
      ValueCell.update(this.blurSecondPassRenderable.values.dOrthographic, orthographic);
    }

    const includeTransparent = SsaoPass.isTransparentEnabled(scene, props);
    if (this.renderable.values.dIncludeTransparent.ref.value !== includeTransparent) {
      needsUpdateSsao = true;

      ValueCell.update(this.renderable.values.dIncludeTransparent, includeTransparent);
    }

    if (this.renderable.values.dIllumination.ref.value !== illuminationMode) {
      needsUpdateSsao = true;

      ValueCell.update(this.renderable.values.dIllumination, illuminationMode);
    }

    if (this.nSamples !== Math.min(props.samples, MaxSamples)) {
      needsUpdateSsao = true;

      this.nSamples = Math.min(props.samples, MaxSamples);
      ValueCell.update(this.renderable.values.uSamples, getSamples(this.nSamples));
      ValueCell.updateIfChanged(this.renderable.values.dNSamples, this.nSamples);
    }

    const multiScale = props.multiScale.name === 'on';
    if (this.renderable.values.dMultiScale.ref.value !== multiScale) {
      needsUpdateSsao = true;

      ValueCell.update(this.renderable.values.dMultiScale, multiScale);
    }

    if (props.multiScale.name === 'on') {
      const mp = props.multiScale.params;
      if (this.levelsCameraScale !== camera.scale || !deepEqual(this.levels, mp.levels)) {
        needsUpdateSsao = true;

        this.levelsCameraScale = camera.scale;
        this.levels = mp.levels;
        const levels = getLevels(mp.levels, camera.scale);
        ValueCell.updateIfChanged(this.renderable.values.dLevels, levels.count);

        ValueCell.update(this.renderable.values.uLevelRadius, levels.radius);
        ValueCell.update(this.renderable.values.uLevelBias, levels.bias);
      }
      ValueCell.updateIfChanged(this.renderable.values.uNearThreshold, mp.nearThreshold);
      ValueCell.updateIfChanged(this.renderable.values.uFarThreshold, mp.farThreshold);
    } else {
      ValueCell.updateIfChanged(this.renderable.values.uRadius, Math.pow(2, props.radius) * camera.scale);
    }
    ValueCell.updateIfChanged(this.renderable.values.uBias, props.bias);

    if (this.blurKernelSize !== Math.min(props.blurKernelSize, MaxBlurKernelSize)) {
      needsUpdateSsaoBlur = true;

      this.blurKernelSize = Math.min(props.blurKernelSize, MaxBlurKernelSize);
      const kernel = getBlurKernel(this.blurKernelSize);

      ValueCell.update(this.blurFirstPassRenderable.values.uKernel, kernel);
      ValueCell.update(this.blurSecondPassRenderable.values.uKernel, kernel);
      ValueCell.update(this.blurFirstPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
      ValueCell.update(this.blurSecondPassRenderable.values.dOcclusionKernelSize, this.blurKernelSize);
    }

    this.resolutionScale = props.resolutionScale;
    const ssaoScale = this.calcSsaoScale(props.resolutionScale);
    if (this.ssaoScale !== ssaoScale) {
      needsUpdateSsao = true;
      needsUpdateDepthHalf = true;

      this.ssaoScale = ssaoScale;

      const sw = Math.floor(w * this.ssaoScale);
      const sh = Math.floor(h * this.ssaoScale);
      this.ssaoDepthTexture.define(sw, sh);
      this.ssaoDepthTransparentTexture.define(sw, sh);
      this.depthBlurProxyTexture.define(sw, sh);

      const hw = Math.floor(sw * 0.5);
      const hh = Math.floor(sh * 0.5);
      const qw = Math.floor(sw * 0.25);
      const qh = Math.floor(sh * 0.25);

      this.depthHalfTargetOpaque.setSize(hw, hh);
      this.depthQuarterTargetOpaque.setSize(qw, qh);

      const depthTexture = this.getDepthTexture();
      ValueCell.update(this.depthHalfRenderableOpaque.values.tColor, depthTexture);
      ValueCell.update(this.renderable.values.tDepth, depthTexture);
      this.updatePackedDepth();

      ValueCell.update(this.renderable.values.tDepthHalf, this.depthHalfTargetOpaque.texture);
      ValueCell.update(this.renderable.values.tDepthQuarter, this.depthQuarterTargetOpaque.texture);

      ValueCell.update(
        this.depthHalfRenderableOpaque.values.uTexSize,
        Vec2.set(this.depthHalfRenderableOpaque.values.uTexSize.ref.value, hw, hh),
      );
      ValueCell.update(
        this.depthQuarterRenderableOpaque.values.uTexSize,
        Vec2.set(this.depthQuarterRenderableOpaque.values.uTexSize.ref.value, qw, qh),
      );

      this.depthHalfTargetTransparent.setSize(hw, hh);
      this.depthQuarterTargetTransparent.setSize(qw, qh);

      const transparentDepthTexture = this.getTransparentDepthTexture();
      ValueCell.update(this.depthHalfRenderableTransparent.values.tColor, transparentDepthTexture);
      ValueCell.update(this.renderable.values.tDepthTransparent, transparentDepthTexture);

      ValueCell.update(this.renderable.values.tDepthHalfTransparent, this.depthHalfTargetTransparent.texture);
      ValueCell.update(this.renderable.values.tDepthQuarterTransparent, this.depthQuarterTargetTransparent.texture);

      ValueCell.update(
        this.depthHalfRenderableTransparent.values.uTexSize,
        Vec2.set(this.depthHalfRenderableTransparent.values.uTexSize.ref.value, hw, hh),
      );
      ValueCell.update(
        this.depthQuarterRenderableTransparent.values.uTexSize,
        Vec2.set(this.depthQuarterRenderableTransparent.values.uTexSize.ref.value, qw, qh),
      );

      ValueCell.update(this.renderable.values.uTexSize, Vec2.set(this.renderable.values.uTexSize.ref.value, sw, sh));
      ValueCell.update(
        this.blurFirstPassRenderable.values.uTexSize,
        Vec2.set(this.blurFirstPassRenderable.values.uTexSize.ref.value, sw, sh),
      );
      ValueCell.update(
        this.blurSecondPassRenderable.values.uTexSize,
        Vec2.set(this.blurSecondPassRenderable.values.uTexSize.ref.value, sw, sh),
      );
    }

    if (needsUpdateSsao) {
      this.renderable.update();
    }

    if (needsUpdateSsaoBlur) {
      this.blurFirstPassRenderable.update();
      this.blurSecondPassRenderable.update();
    }

    if (needsUpdateDepthHalf) {
      this.depthHalfRenderableOpaque.update();
      this.depthHalfRenderableTransparent.update();
    }
  }

  render(b: PassBuilder, camera: ICamera) {
    const { graphics } = this;
    const timing = isTimingMode && b.drawing;
    if (timing) graphics.timer.group('SSAO.render');

    const { x, y, width, height } = camera.viewport;
    const state = RenderState.Default;

    const includeTransparent = this.renderable.values.dIncludeTransparent.ref.value;
    const multiScale = this.renderable.values.dMultiScale.ref.value;

    const sx = Math.floor(x * this.ssaoScale);
    const sy = Math.floor(y * this.ssaoScale);
    const sw = Math.ceil(width * this.ssaoScale);
    const sh = Math.ceil(height * this.ssaoScale);

    if (b.drawing) {
      graphics.setViewport(sx, sy, sw, sh);
      graphics.setScissor(sx, sy, sw, sh);
    }

    if (multiScale) {
      // half-resolution viewport (matches dimensions of depthHalfTarget*)
      if (b.drawing) {
        const hsx = Math.floor(sx * 0.5);
        const hsy = Math.floor(sy * 0.5);
        const hsw = Math.ceil(sw * 0.5);
        const hsh = Math.ceil(sh * 0.5);
        graphics.setViewport(hsx, hsy, hsw, hsh);
        graphics.setScissor(hsx, hsy, hsw, hsh);
      }

      if (timing) graphics.timer.group('SSAO.half');
      b.begin(this.depthHalfTargetOpaque, { label: 'ssao half' });
      b.quad(this.depthHalfRenderableOpaque, state);
      if (includeTransparent) {
        b.begin(this.depthHalfTargetTransparent, { label: 'ssao half transparent' });
        b.quad(this.depthHalfRenderableTransparent, state);
      }
      if (timing) graphics.timer.groupEnd('SSAO.half');

      // quarter-resolution viewport (matches dimensions of depthQuarterTarget*)
      if (b.drawing) {
        const qsx = Math.floor(sx * 0.25);
        const qsy = Math.floor(sy * 0.25);
        const qsw = Math.ceil(sw * 0.25);
        const qsh = Math.ceil(sh * 0.25);
        graphics.setViewport(qsx, qsy, qsw, qsh);
        graphics.setScissor(qsx, qsy, qsw, qsh);
      }

      if (timing) graphics.timer.group('SSAO.quarter');
      b.begin(this.depthQuarterTargetOpaque, { label: 'ssao quarter' });
      b.quad(this.depthQuarterRenderableOpaque, state);
      if (includeTransparent) {
        b.begin(this.depthQuarterTargetTransparent, { label: 'ssao quarter transparent' });
        b.quad(this.depthQuarterRenderableTransparent, state);
      }
      if (timing) graphics.timer.groupEnd('SSAO.quarter');

      // restore full-scale viewport for SSAO + blur passes
      if (b.drawing) {
        graphics.setViewport(sx, sy, sw, sh);
        graphics.setScissor(sx, sy, sw, sh);
      }
    }

    this.ssaoDepthTexture.attachFramebuffer(this.framebuffer, 'color0');
    ValueCell.update(this.renderable.values.uTransparencyFlag, 0);
    b.begin(this.framebuffer, { label: 'ssao' });
    b.quad(this.renderable, state);

    if (this.blur) {
      if (timing) graphics.timer.group('SSAO.blurOpaque');
      ValueCell.update(this.blurFirstPassRenderable.values.tSsaoDepth, this.ssaoDepthTexture);
      this.blurFirstPassRenderable.update();
      b.begin(this.blurFirstPassFramebuffer, { label: 'ssao blur' });
      b.quad(this.blurFirstPassRenderable, state);

      this.ssaoDepthTexture.attachFramebuffer(this.blurSecondPassFramebuffer, 'color0');
      b.begin(this.blurSecondPassFramebuffer, { label: 'ssao blur' });
      b.quad(this.blurSecondPassRenderable, state);
      if (timing) graphics.timer.groupEnd('SSAO.blurOpaque');
    }
    if (includeTransparent) {
      this.ssaoDepthTransparentTexture.attachFramebuffer(this.framebuffer, 'color0');
      ValueCell.update(this.renderable.values.uTransparencyFlag, 1);
      b.begin(this.framebuffer, { label: 'ssao transparent' });
      b.quad(this.renderable, state);

      if (this.blur) {
        if (timing) graphics.timer.group('SSAO.blurTransparent ');
        ValueCell.update(this.blurFirstPassRenderable.values.tSsaoDepth, this.ssaoDepthTransparentTexture);
        this.blurFirstPassRenderable.update();
        b.begin(this.blurFirstPassFramebuffer, { label: 'ssao blur transparent' });
        b.quad(this.blurFirstPassRenderable, state);

        this.ssaoDepthTransparentTexture.attachFramebuffer(this.blurSecondPassFramebuffer, 'color0');
        b.begin(this.blurSecondPassFramebuffer, { label: 'ssao blur transparent' });
        b.quad(this.blurSecondPassRenderable, state);
        if (timing) graphics.timer.groupEnd('SSAO.blurTransparent ');
      }
    }
    if (timing) graphics.timer.groupEnd('SSAO.render');
  }
}

export const SsaoSchema = {
  ...QuadSchema,
  tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tDepthHalf: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tDepthQuarter: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  dPackedDepth: DefineSpec('boolean'),
  dPackedDepthLevels: DefineSpec('boolean'),

  dIllumination: DefineSpec('boolean'),
  uTransparencyFlag: UniformSpec('i'),
  dIncludeTransparent: DefineSpec('boolean'),
  tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tDepthHalfTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tDepthQuarterTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),

  uSamples: UniformSpec('v3[]', undefined, undefined, MaxSamples, 'dNSamples'),
  dNSamples: DefineSpec('number'),

  uProjection: { ...UniformSpec('m4'), depthRange: 'projection' as const },
  uInvProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },
  uBounds: { ...UniformSpec('v4'), windowRect: 'bounds' as const },
  /** the near and far plane, for the linear depth of the output */
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),

  uTexSize: UniformSpec('v2'),

  uRadius: UniformSpec('f'),
  uBias: UniformSpec('f'),

  dMultiScale: DefineSpec('boolean'),
  dLevels: DefineSpec('number'),
  uLevelRadius: UniformSpec('f[]', undefined, undefined, MaxLevels, 'dLevels'),
  uLevelBias: UniformSpec('f[]', undefined, undefined, MaxLevels, 'dLevels'),
  uNearThreshold: UniformSpec('f'),
  uFarThreshold: UniformSpec('f'),
};
export const SsaoShaderCode = ShaderCode('ssao');

type SsaoRenderable = ComputeRenderable<Values<typeof SsaoSchema>>;

function getSsaoRenderable(
  ctx: GraphicsContext,
  depthTexture: Texture,
  depthHalfTexture: Texture,
  depthQuarterTexture: Texture,
  transparentDepthTexture: Texture,
  transparentDepthHalfTexture: Texture,
  transparentDepthQuarterTexture: Texture,
  packedDepthLevels: boolean,
): SsaoRenderable {
  const drawingBufferSize = ctx.getDrawingBufferSize();
  const values: Values<typeof SsaoSchema> = {
    ...QuadValues,
    tDepth: ValueCell.create(depthTexture),
    tDepthHalf: ValueCell.create(depthHalfTexture),
    tDepthQuarter: ValueCell.create(depthQuarterTexture),
    dPackedDepth: ValueCell.create(false),
    dPackedDepthLevels: ValueCell.create(packedDepthLevels),

    dIllumination: ValueCell.create(false),
    dIncludeTransparent: ValueCell.create(true),
    uTransparencyFlag: ValueCell.create(0),
    tDepthTransparent: ValueCell.create(transparentDepthTexture),
    tDepthHalfTransparent: ValueCell.create(transparentDepthHalfTexture),
    tDepthQuarterTransparent: ValueCell.create(transparentDepthQuarterTexture),

    uSamples: ValueCell.create(getSamples(32)),
    dNSamples: ValueCell.create(32),

    uProjection: ValueCell.create(Mat4.identity()),
    uInvProjection: ValueCell.create(Mat4.identity()),
    uBounds: ValueCell.create(Vec4()),

    uTexSize: ValueCell.create(Vec2.create(drawingBufferSize.width, drawingBufferSize.height)),

    uRadius: ValueCell.create(Math.pow(2, 5)),
    uBias: ValueCell.create(0.8),

    dMultiScale: ValueCell.create(false),
    dLevels: ValueCell.create(3),
    uLevelRadius: ValueCell.create([Math.pow(2, 2), Math.pow(2, 5), Math.pow(2, 8)]),
    uLevelBias: ValueCell.create([0.8, 0.8, 0.8]),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    uNearThreshold: ValueCell.create(10.0),
    uFarThreshold: ValueCell.create(1500.0),
  };

  const schema = { ...SsaoSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', SsaoShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

export const SsaoBlurSchema = {
  ...QuadSchema,
  tSsaoDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),

  uKernel: UniformSpec('f[]', undefined, undefined, MaxBlurKernelSize, 'dOcclusionKernelSize'),
  dOcclusionKernelSize: DefineSpec('number'),
  uBlurDepthBias: UniformSpec('f'),

  uBlurDirectionX: UniformSpec('f'),
  uBlurDirectionY: UniformSpec('f'),

  uInvProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  uBounds: { ...UniformSpec('v4'), windowRect: 'bounds' as const },
  dOrthographic: DefineSpec('number'),
};
export const SsaoBlurShaderCode = ShaderCode('ssao_blur');

type SsaoBlurRenderable = ComputeRenderable<Values<typeof SsaoBlurSchema>>;

function getSsaoBlurRenderable(
  ctx: GraphicsContext,
  ssaoDepthTexture: Texture,
  direction: 'horizontal' | 'vertical',
): SsaoBlurRenderable {
  const values: Values<typeof SsaoBlurSchema> = {
    ...QuadValues,
    tSsaoDepth: ValueCell.create(ssaoDepthTexture),
    uTexSize: ValueCell.create(Vec2.create(ssaoDepthTexture.getWidth(), ssaoDepthTexture.getHeight())),

    uKernel: ValueCell.create(getBlurKernel(15)),
    dOcclusionKernelSize: ValueCell.create(15),
    uBlurDepthBias: ValueCell.create(0.5),

    uBlurDirectionX: ValueCell.create(direction === 'horizontal' ? 1 : 0),
    uBlurDirectionY: ValueCell.create(direction === 'vertical' ? 1 : 0),

    uInvProjection: ValueCell.create(Mat4.identity()),
    uNear: ValueCell.create(0.0),
    uFar: ValueCell.create(10000.0),
    uBounds: ValueCell.create(Vec4()),
    dOrthographic: ValueCell.create(0),
  };

  const schema = { ...SsaoBlurSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', SsaoBlurShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

function getBlurKernel(kernelSize: number): number[] {
  const sigma = kernelSize / 3.0;
  const halfKernelSize = Math.floor((kernelSize + 1) / 2);

  const kernel = [];
  for (let x = 0; x < halfKernelSize; x++) {
    kernel.push((1.0 / (Math.sqrt(2 * Math.PI) * sigma)) * Math.exp((-x * x) / (2 * sigma * sigma)));
  }

  return kernel;
}

const pcg = new PCG();
function getRandomHemisphereVector(): Vec3 {
  const v = Vec3();
  while (true) {
    const x = pcg.float() * 2 - 1;
    const y = pcg.float() * 2 - 1;
    if (x * x + y * y < 1) {
      const z = 2 * Math.sqrt(1 - x * x - y * y) * (pcg.float() < 0.5 ? -1 : 1);
      Vec3.set(v, x, y, z);
      Vec3.normalize(v, v);
      Vec3.scale(v, v, pcg.float());
      break;
    }
  }
  if (v[2] < 0) v[2] = -v[2];
  return v;
}

function generateBlueNoiseVectors(count: number, out: Vec3[]) {
  if (out.length >= count) return out;
  if (out.length === 0) out.push(getRandomHemisphereVector());

  const candidateCount = Math.max(10, Math.min(30, Math.floor(count / 10)));

  for (let i = out.length; i < count; i++) {
    let bestCandidate: Vec3;
    let bestDistance = -1;

    for (let j = 0; j < candidateCount; j++) {
      const candidate = getRandomHemisphereVector();

      let minDistance = Infinity;
      for (const existingVector of out) {
        const distance = Vec3.distance(candidate, existingVector);
        minDistance = Math.min(minDistance, distance);
      }

      if (minDistance > bestDistance) {
        bestDistance = minDistance;
        bestCandidate = candidate;
      }
    }

    out.push(bestCandidate!);
  }

  return out;
}

let _RandomHemisphereVectors: Vec3[] = [];
function getRandomHemisphereVectors(count: number): Vec3[] {
  if (_RandomHemisphereVectors.length < count) {
    _RandomHemisphereVectors = generateBlueNoiseVectors(count, _RandomHemisphereVectors);
  }
  return _RandomHemisphereVectors;
}

function getSamples(nSamples: number): number[] {
  const rhv = getRandomHemisphereVectors(nSamples);
  const samples = [];
  for (let i = 0; i < nSamples; i++) {
    let scale = (i * i + 2.0 * i + 1) / (nSamples * nSamples);
    scale = 0.1 + scale * (1.0 - 0.1);

    samples.push(rhv[i][0] * scale);
    samples.push(rhv[i][1] * scale);
    samples.push(rhv[i][2] * scale);
  }

  return samples;
}
