/**
 * Copyright (c) 2024-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 *
 * Partially adapted from three.js, The MIT License, Copyright © 2010-2024 three.js authors
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { type ComputeRenderable, createComputeRenderable } from '@molstar/graphics/gl/renderable';
import { TextureSpec, UniformSpec, DefineSpec, type Values } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Scene } from '@molstar/graphics/gl/scene';
import { type Texture, createNullTexture } from '@molstar/graphics/gl/texture';
import { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Viewport } from '../camera/util.js';
import type { ICamera } from '../camera.js';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import { isTimingMode } from '@molstar/core/util/debug';
import { memoize1 } from '@molstar/core/util/memoize';
import type { PostprocessingProps } from './postprocessing.js';
import type { PassTargets } from './pass-targets.js';

const MipCount = 5;

export const BloomParams = {
  strength: PD.Numeric(1, { min: 0, max: 3, step: 0.1 }),
  radius: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
  threshold: PD.Numeric(
    0,
    { min: 0, max: 1, step: 0.01 },
    { description: 'Luminosity threshold', hideIf: (p) => p.mode === 'emissive' },
  ),
  mode: PD.Select('emissive', [
    ['luminosity', 'Luminosity'],
    ['emissive', 'Emissive'],
  ] as const),
  transparency: PD.Boolean(true, {
    description: 'Include transparent objects in emissive bloom',
    hideIf: (p) => p.mode !== 'emissive',
  }),
};
export type BloomProps = PD.Values<typeof BloomParams>;

export class BloomPass {
  private readonly graphics: GraphicsContext;
  static isEnabled(props: PostprocessingProps) {
    return props.enabled && props.bloom.name === 'on';
  }

  /** Whether bloom draws: enabled, and in emissive mode (the default) something in the scene is emissive */
  static isActive(scene: Scene, props: PostprocessingProps) {
    return (
      props.enabled &&
      props.bloom.name === 'on' &&
      !(props.bloom.params.mode === 'emissive' && scene.emissiveAverage === 0)
    );
  }

  readonly emissiveTarget: RenderTarget;
  readonly compositeTarget: RenderTarget;

  private readonly luminosityTarget: RenderTarget;
  private readonly horizontalBlurTargets: RenderTarget[] = [];
  private readonly verticalBlurTargets: RenderTarget[] = [];

  private readonly luminosityRenderable: LuminosityRenderable;
  private readonly blurRenderable: BlurRenderable;
  private readonly compositeRenderable: CompositeRenderable;

  constructor(targets: PassTargets, width: number, height: number) {
    const graphics = (this.graphics = targets.graphics);
    this.emissiveTarget = targets.createRenderTarget(width, height, 'depth-stencil', 'uint8', 'linear');

    this.luminosityTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');
    this.compositeTarget = targets.createRenderTarget(width, height, 'none', 'uint8', 'linear');

    let blurWidth = Math.round(width / 2);
    let blurHeight = Math.round(height / 2);
    for (let i = 0; i < MipCount; ++i) {
      this.horizontalBlurTargets[i] = targets.createRenderTarget(blurWidth, blurHeight, 'none', 'uint8', 'linear');
      this.verticalBlurTargets[i] = targets.createRenderTarget(blurWidth, blurHeight, 'none', 'uint8', 'linear');
      blurWidth = Math.round(blurWidth / 2);
      blurHeight = Math.round(blurHeight / 2);
    }

    const nullTexture = createNullTexture();
    this.luminosityRenderable = getLuminosityRenderable(
      graphics,
      nullTexture,
      nullTexture,
      nullTexture,
      nullTexture,
      nullTexture,
    );
    this.blurRenderable = getBlurRenderable(graphics, nullTexture);
    this.compositeRenderable = getCompositeRenderable(
      graphics,
      width,
      height,
      this.verticalBlurTargets[0].texture,
      this.verticalBlurTargets[1].texture,
      this.verticalBlurTargets[2].texture,
      this.verticalBlurTargets[3].texture,
      this.verticalBlurTargets[4].texture,
    );
  }

  getByteCount() {
    return (
      this.emissiveTarget.getByteCount() +
      this.luminosityTarget.getByteCount() +
      this.compositeTarget.getByteCount() +
      this.horizontalBlurTargets.reduce((sum, t) => sum + t.getByteCount(), 0) +
      this.verticalBlurTargets.reduce((sum, t) => sum + t.getByteCount(), 0)
    );
  }

  setSize(width: number, height: number) {
    const w = this.luminosityTarget.getWidth();
    const h = this.luminosityTarget.getHeight();

    if (width !== w || height !== h) {
      this.emissiveTarget.setSize(width, height);
      this.luminosityTarget.setSize(width, height);
      this.compositeTarget.setSize(width, height);

      let blurWidth = Math.round(width / 2);
      let blurHeight = Math.round(height / 2);
      for (let i = 0; i < MipCount; ++i) {
        this.horizontalBlurTargets[i].setSize(blurWidth, blurHeight);
        this.verticalBlurTargets[i].setSize(blurWidth, blurHeight);
        blurWidth = Math.round(blurWidth / 2);
        blurHeight = Math.round(blurHeight / 2);
      }

      ValueCell.update(
        this.luminosityRenderable.values.uTexSizeInv,
        Vec2.set(this.luminosityRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height),
      );
      ValueCell.update(
        this.compositeRenderable.values.uTexSizeInv,
        Vec2.set(this.compositeRenderable.values.uTexSizeInv.ref.value, 1 / width, 1 / height),
      );
    }
  }

  update(
    colorOpaque: Texture,
    colorTransparent: Texture,
    emissive: Texture,
    depthOpaque: Texture,
    depthTransparent: Texture,
    props: BloomProps,
    camera: ICamera,
    opaqueFogged: boolean,
    /** whether the transparent color and depth are composited in (read), as in the postprocessing compose */
    blendTransparency: boolean,
  ) {
    let luminosityNeedsUpdate = false;

    if (this.luminosityRenderable.values.tColorOpaque.ref.value !== colorOpaque) {
      ValueCell.update(this.luminosityRenderable.values.tColorOpaque, colorOpaque);
      luminosityNeedsUpdate = true;
    }

    if (this.luminosityRenderable.values.tColorTransparent.ref.value !== colorTransparent) {
      ValueCell.update(this.luminosityRenderable.values.tColorTransparent, colorTransparent);
      luminosityNeedsUpdate = true;
    }

    if (this.luminosityRenderable.values.tEmissive.ref.value !== emissive) {
      ValueCell.update(this.luminosityRenderable.values.tEmissive, emissive);
      luminosityNeedsUpdate = true;
    }

    if (this.luminosityRenderable.values.tDepthOpaque.ref.value !== depthOpaque) {
      ValueCell.update(this.luminosityRenderable.values.tDepthOpaque, depthOpaque);
      luminosityNeedsUpdate = true;
    }

    if (this.luminosityRenderable.values.tDepthTransparent.ref.value !== depthTransparent) {
      ValueCell.update(this.luminosityRenderable.values.tDepthTransparent, depthTransparent);
      luminosityNeedsUpdate = true;
    }

    if (this.luminosityRenderable.values.dBlendTransparency.ref.value !== blendTransparency) {
      ValueCell.update(this.luminosityRenderable.values.dBlendTransparency, blendTransparency);
      luminosityNeedsUpdate = true;
    }

    if (this.luminosityRenderable.values.dMode.ref.value !== props.mode) {
      ValueCell.update(this.luminosityRenderable.values.dMode, props.mode);
      luminosityNeedsUpdate = true;
    }

    ValueCell.updateIfChanged(this.luminosityRenderable.values.uLuminosityThreshold, props.threshold);
    ValueCell.updateIfChanged(this.luminosityRenderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.luminosityRenderable.values.uFar, camera.far);
    ValueCell.updateIfChanged(this.luminosityRenderable.values.uIsOrtho, camera.state.mode === 'orthographic' ? 1 : 0);
    ValueCell.updateIfChanged(this.luminosityRenderable.values.uFogNear, camera.fogNear);
    ValueCell.updateIfChanged(this.luminosityRenderable.values.uFogFar, camera.fogFar);
    ValueCell.updateIfChanged(this.luminosityRenderable.values.uOpaqueFogged, opaqueFogged);

    if (luminosityNeedsUpdate) {
      this.luminosityRenderable.update();
    }

    ValueCell.updateIfChanged(this.compositeRenderable.values.uBloomRadius, props.radius);
    ValueCell.updateIfChanged(this.compositeRenderable.values.uBloomStrength, props.strength);
  }

  /** Writes bloom to `compositeTarget`; caller blends it. */
  render(b: PassBuilder, viewport: Viewport) {
    const { graphics } = this;
    const timing = isTimingMode && b.drawing;
    if (timing) graphics.timer.group('BloomPass.render');
    const { x, y, width, height } = viewport;
    const state = RenderState.Default;
    if (b.drawing) {
      graphics.setViewport(x, y, width, height);
      graphics.setScissor(x, y, width, height);
    }

    b.begin(this.luminosityTarget, { label: 'bloom luminosity' });
    b.quad(this.luminosityRenderable, state);

    for (let i = 0; i < MipCount; ++i) {
      const blurWidth = this.horizontalBlurTargets[i].getWidth();
      const blurHeight = this.horizontalBlurTargets[i].getHeight();
      if (b.drawing) {
        graphics.setViewport(0, 0, blurWidth, blurHeight);
        graphics.setScissor(0, 0, blurWidth, blurHeight);
      }

      ValueCell.update(this.blurRenderable.values.dKernelRadius, BlurKernelSizes[i]);
      ValueCell.update(this.blurRenderable.values.uGaussianCoefficients, getBlurCoefficients(BlurKernelSizes[i]));
      ValueCell.update(
        this.blurRenderable.values.uTexSizeInv,
        Vec2.set(this.blurRenderable.values.uTexSizeInv.ref.value, 1 / blurWidth, 1 / blurHeight),
      );

      b.begin(this.horizontalBlurTargets[i], { label: 'bloom blur' });
      ValueCell.update(
        this.blurRenderable.values.tInput,
        i === 0 ? this.luminosityTarget.texture : this.verticalBlurTargets[i - 1].texture,
      );
      ValueCell.update(this.blurRenderable.values.uDirection, BlurDirectionX);
      this.blurRenderable.update();
      b.quad(this.blurRenderable, state);

      b.begin(this.verticalBlurTargets[i], { label: 'bloom blur' });
      ValueCell.update(this.blurRenderable.values.tInput, this.horizontalBlurTargets[i].texture);
      ValueCell.update(this.blurRenderable.values.uDirection, BlurDirectionY);
      this.blurRenderable.update();
      b.quad(this.blurRenderable, state);
    }

    if (b.drawing) {
      graphics.setViewport(x, y, width, height);
      graphics.setScissor(x, y, width, height);
    }

    b.begin(this.compositeTarget, { label: 'bloom composite' });
    this.compositeRenderable.update();
    b.quad(this.compositeRenderable, state);

    if (timing) graphics.timer.groupEnd('BloomPass.render');
  }
}

//

export const LuminositySchema = {
  ...QuadSchema,
  tColorOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tColorTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tEmissive: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSizeInv: UniformSpec('v2'),

  uDefaultColor: UniformSpec('v3'),
  uDefaultOpacity: UniformSpec('f'),
  uLuminosityThreshold: UniformSpec('f'),
  uSmoothWidth: UniformSpec('f'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  uIsOrtho: UniformSpec('f'),
  uFogNear: UniformSpec('f'),
  uFogFar: UniformSpec('f'),
  uOpaqueFogged: UniformSpec('b'),
  dBlendTransparency: DefineSpec('boolean'),
  dMode: DefineSpec('string', ['luminosity', 'emissive']),
};
export const LuminosityShaderCode = ShaderCode('Bloom Luminosity');
type LuminosityRenderable = ComputeRenderable<Values<typeof LuminositySchema>>;

function getLuminosityRenderable(
  ctx: GraphicsContext,
  colorOpaqueTexture: Texture,
  colorTransparentTexture: Texture,
  emissiveTexture: Texture,
  depthOpaqueTexture: Texture,
  depthTransparentTexture: Texture,
): LuminosityRenderable {
  const width = colorOpaqueTexture.getWidth();
  const height = colorOpaqueTexture.getHeight();

  const values: Values<typeof LuminositySchema> = {
    ...QuadValues,
    tColorOpaque: ValueCell.create(colorOpaqueTexture),
    tColorTransparent: ValueCell.create(colorTransparentTexture),
    tEmissive: ValueCell.create(emissiveTexture),
    tDepthOpaque: ValueCell.create(depthOpaqueTexture),
    tDepthTransparent: ValueCell.create(depthTransparentTexture),
    uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

    uDefaultColor: ValueCell.create(Vec3()),
    uDefaultOpacity: ValueCell.create(0),
    uLuminosityThreshold: ValueCell.create(0),
    uSmoothWidth: ValueCell.create(1),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    uIsOrtho: ValueCell.create(0),
    uFogNear: ValueCell.create(10000),
    uFogFar: ValueCell.create(10000),
    uOpaqueFogged: ValueCell.create(true),
    dBlendTransparency: ValueCell.create(true),
    dMode: ValueCell.create('emissive'),
  };

  const schema = { ...LuminositySchema };
  const renderItem = ctx.createComputeRenderItem('triangles', LuminosityShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

function _getBlurCoefficients(kernelRadius: number) {
  const coefficients: number[] = [];
  for (let i = 0; i < kernelRadius; ++i) {
    coefficients.push((0.39894 * Math.exp((-0.5 * i * i) / (kernelRadius * kernelRadius))) / kernelRadius);
  }
  return coefficients;
}
const getBlurCoefficients = memoize1(_getBlurCoefficients);

const BlurKernelSizes = [3, 5, 7, 9, 11];

const BlurDirectionX = Vec2.create(1, 0);
const BlurDirectionY = Vec2.create(0, 1);

export const BlurSchema = {
  ...QuadSchema,
  tInput: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  uTexSizeInv: UniformSpec('v2'),

  uDirection: UniformSpec('v2'),
  uGaussianCoefficients: UniformSpec('f[]', undefined, undefined, Math.max(...BlurKernelSizes), 'dKernelRadius'),
  dKernelRadius: DefineSpec('number'),
};
export const BlurShaderCode = ShaderCode('Bloom Blur');
type BlurRenderable = ComputeRenderable<Values<typeof BlurSchema>>;

function getBlurRenderable(ctx: GraphicsContext, inputTexture: Texture): BlurRenderable {
  const width = inputTexture.getWidth();
  const height = inputTexture.getHeight();

  const values: Values<typeof BlurSchema> = {
    ...QuadValues,
    tInput: ValueCell.create(inputTexture),
    uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

    uDirection: ValueCell.create(Vec2()),
    uGaussianCoefficients: ValueCell.create([]),
    dKernelRadius: ValueCell.create(BlurKernelSizes[0]),
  };

  const schema = { ...BlurSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', BlurShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

export const CompositeSchema = {
  ...QuadSchema,
  tBlur1: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tBlur2: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tBlur3: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tBlur4: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  tBlur5: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
  uTexSizeInv: UniformSpec('v2'),

  uBloomStrength: UniformSpec('f'),
  uBloomRadius: UniformSpec('f'),
  uBloomFactors: UniformSpec('f[]', undefined, undefined, 5),
  uBloomTints: UniformSpec('v3[]', undefined, undefined, 5),
};
export const CompositeShaderCode = ShaderCode('Bloom Composite');
type CompositeRenderable = ComputeRenderable<Values<typeof CompositeSchema>>;

function getCompositeRenderable(
  ctx: GraphicsContext,
  width: number,
  height: number,
  blurTexture1: Texture,
  blurTexture2: Texture,
  blurTexture3: Texture,
  blurTexture4: Texture,
  blurTexture5: Texture,
): CompositeRenderable {
  const values: Values<typeof CompositeSchema> = {
    ...QuadValues,
    uTexSizeInv: ValueCell.create(Vec2.create(1 / width, 1 / height)),

    tBlur1: ValueCell.create(blurTexture1),
    tBlur2: ValueCell.create(blurTexture2),
    tBlur3: ValueCell.create(blurTexture3),
    tBlur4: ValueCell.create(blurTexture4),
    tBlur5: ValueCell.create(blurTexture5),

    uBloomStrength: ValueCell.create(1),
    uBloomRadius: ValueCell.create(0),
    uBloomFactors: ValueCell.create([1.0, 0.8, 0.6, 0.4, 0.2]),
    uBloomTints: ValueCell.create(new Array(5 * 3).fill(1)),
  };

  const schema = { ...CompositeSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', CompositeShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
