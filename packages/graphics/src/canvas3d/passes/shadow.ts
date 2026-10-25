/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Áron Samuel Kovács <aron.kovacs@mail.muni.cz>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { PassBuilder } from './pass-builder.js';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import {
  TextureSpec,
  type Values,
  UniformSpec,
  DefineSpec,
  MaxLightCount,
} from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import { ValueCell } from '@molstar/core/util';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { ICamera } from '../camera.js';
import { getTransformedLightDirection, type Light } from '@molstar/graphics/gl/renderer';
import type { PostprocessingProps } from './postprocessing.js';
import type { PassTargets } from './pass-targets.js';

export const ShadowParams = {
  steps: PD.Numeric(1, { min: 1, max: 64, step: 1 }),
  maxDistance: PD.Numeric(3, { min: 0, max: 256, step: 1 }),
  tolerance: PD.Numeric(1.0, { min: 0.0, max: 10.0, step: 0.1 }),
};

export type ShadowProps = PD.Values<typeof ShadowParams>;

export class ShadowPass {
  readonly graphics: GraphicsContext;
  static isEnabled(props: PostprocessingProps) {
    return props.enabled && props.shadow.name !== 'off';
  }

  readonly target: RenderTarget;
  private readonly renderable: ShadowsRenderable;

  private invProjection = Mat4.identity();
  private invHeadRotation = Mat4.identity();

  constructor(targets: PassTargets, width: number, height: number, depthTextureOpaque: Texture) {
    const graphics = (this.graphics = targets.graphics);
    this.target = targets.createRenderTarget(width, height, 'none');
    this.renderable = getShadowsRenderable(graphics, depthTextureOpaque);
  }

  getByteCount() {
    return this.target.getByteCount();
  }

  setSize(width: number, height: number) {
    const [w, h] = this.renderable.values.uTexSize.ref.value;
    if (width !== w || height !== h) {
      this.target.setSize(width, height);
      ValueCell.update(
        this.renderable.values.uTexSize,
        Vec2.set(this.renderable.values.uTexSize.ref.value, width, height),
      );
    }
  }

  update(camera: ICamera, light: Light, ambientColor: Vec3, props: ShadowProps) {
    let needsUpdateShadows = false;

    const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

    const [w, h] = this.renderable.values.uTexSize.ref.value;
    const v = camera.viewport;

    ValueCell.update(this.renderable.values.uProjection, camera.projection);
    ValueCell.update(this.renderable.values.uInvProjection, Mat4.invert(this.invProjection, camera.projection));

    Vec4.set(this.renderable.values.uBounds.ref.value, v.x / w, v.y / h, (v.x + v.width) / w, (v.y + v.height) / h);
    ValueCell.update(this.renderable.values.uBounds, this.renderable.values.uBounds.ref.value);

    ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
    if (this.renderable.values.dOrthographic.ref.value !== orthographic) {
      ValueCell.update(this.renderable.values.dOrthographic, orthographic);
      needsUpdateShadows = true;
    }

    ValueCell.updateIfChanged(this.renderable.values.uMaxDistance, props.maxDistance * camera.scale);
    ValueCell.updateIfChanged(this.renderable.values.uTolerance, props.tolerance * camera.scale);
    if (this.renderable.values.dSteps.ref.value !== props.steps) {
      ValueCell.update(this.renderable.values.dSteps, props.steps);
      needsUpdateShadows = true;
    }

    const hasHeadRotation = !Mat4.isZero(camera.headRotation);
    if (hasHeadRotation) {
      ValueCell.update(
        this.renderable.values.uLightDirection,
        getTransformedLightDirection(light, Mat4.invert(this.invHeadRotation, camera.headRotation)),
      );
    } else {
      ValueCell.update(this.renderable.values.uLightDirection, light.direction);
    }
    ValueCell.update(this.renderable.values.uLightColor, light.color);
    if (this.renderable.values.dLightCount.ref.value !== light.count) {
      ValueCell.update(this.renderable.values.dLightCount, light.count);
      needsUpdateShadows = true;
    }
    ValueCell.update(this.renderable.values.uAmbientColor, ambientColor);

    if (needsUpdateShadows) {
      this.renderable.update();
    }
  }

  render(b: PassBuilder) {
    b.begin(this.target, { label: 'shadows' });
    b.quad(this.renderable, RenderState.Default);
  }
}

export const ShadowsSchema = {
  ...QuadSchema,
  tDepth: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),

  uProjection: { ...UniformSpec('m4'), depthRange: 'projection' as const },
  uInvProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },
  uBounds: { ...UniformSpec('v4'), windowRect: 'bounds' as const },

  dOrthographic: DefineSpec('number'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),

  dSteps: DefineSpec('number'),
  uMaxDistance: UniformSpec('f'),
  uTolerance: UniformSpec('f'),

  uLightDirection: UniformSpec('v3[]', undefined, undefined, MaxLightCount, 'dLightCount'),
  uLightColor: UniformSpec('v3[]', undefined, undefined, MaxLightCount, 'dLightCount'),
  dLightCount: DefineSpec('number'),
  uAmbientColor: UniformSpec('v3'),
};
export const ShadowsShaderCode = ShaderCode('shadows');
type ShadowsRenderable = ComputeRenderable<Values<typeof ShadowsSchema>>;

function getShadowsRenderable(ctx: GraphicsContext, depthTexture: Texture): ShadowsRenderable {
  const width = depthTexture.getWidth();
  const height = depthTexture.getHeight();

  const values: Values<typeof ShadowsSchema> = {
    ...QuadValues,
    tDepth: ValueCell.create(depthTexture),
    uTexSize: ValueCell.create(Vec2.create(width, height)),

    uProjection: ValueCell.create(Mat4.identity()),
    uInvProjection: ValueCell.create(Mat4.identity()),
    uBounds: ValueCell.create(Vec4()),

    dOrthographic: ValueCell.create(0),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),

    dSteps: ValueCell.create(1),
    uMaxDistance: ValueCell.create(3.0),
    uTolerance: ValueCell.create(1.0),

    uLightDirection: ValueCell.create([]),
    uLightColor: ValueCell.create([]),
    dLightCount: ValueCell.create(0),
    uAmbientColor: ValueCell.create(Vec3()),
  };

  const schema = { ...ShadowsSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', ShadowsShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
