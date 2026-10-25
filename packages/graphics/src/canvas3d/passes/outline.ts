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
import { TextureSpec, type Values, UniformSpec, DefineSpec } from '@molstar/graphics/gl/renderable/schema';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture } from '@molstar/graphics/gl/texture';
import { ValueCell } from '@molstar/core/util';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import { Mat4, Vec2 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { RenderTarget } from '@molstar/graphics/gl/render-target';
import type { ICamera } from '../camera.js';
import { Color } from '@molstar/core/util/color';
import type { PostprocessingProps } from './postprocessing.js';
import type { PassTargets } from './pass-targets.js';

export const OutlineParams = {
  scale: PD.Numeric(1, { min: 1, max: 5, step: 1 }),
  threshold: PD.Numeric(0.33, { min: 0.01, max: 1, step: 0.01 }),
  color: PD.Color(Color(0x000000)),
  includeTransparent: PD.Boolean(true, { description: 'Whether to show outline for transparent objects' }),
};

export type OutlineProps = PD.Values<typeof OutlineParams>;

export class OutlinePass {
  private readonly graphics: GraphicsContext;
  static isEnabled(props: PostprocessingProps) {
    return props.enabled && props.outline.name !== 'off';
  }

  readonly target: RenderTarget;
  private readonly renderable: OutlinesRenderable;

  constructor(
    targets: PassTargets,
    width: number,
    height: number,
    depthTextureTransparent: Texture,
    depthTextureOpaque: Texture,
  ) {
    const graphics = (this.graphics = targets.graphics);
    this.target = targets.createRenderTarget(width, height, 'none');
    this.renderable = getOutlinesRenderable(graphics, depthTextureOpaque, depthTextureTransparent, true);
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

  update(camera: ICamera, props: OutlineProps, depthTextureTransparent: Texture, depthTextureOpaque: Texture) {
    let needsUpdate = false;

    const orthographic = camera.state.mode === 'orthographic' ? 1 : 0;

    const invProjection = this.renderable.values.uInvProjection.ref.value;
    Mat4.invert(invProjection, camera.projection);

    const transparentOutline = props.includeTransparent ?? true;
    const outlineScale = Math.max(1, Math.round(props.scale * this.graphics.pixelRatio)) - 1;
    const outlineThreshold = 50 * props.threshold * this.graphics.pixelRatio;

    ValueCell.updateIfChanged(this.renderable.values.uNear, camera.near);
    ValueCell.updateIfChanged(this.renderable.values.uFar, camera.far);
    ValueCell.update(this.renderable.values.uInvProjection, invProjection);
    if (this.renderable.values.dTransparentOutline.ref.value !== transparentOutline) {
      needsUpdate = true;
      ValueCell.update(this.renderable.values.dTransparentOutline, transparentOutline);
    }
    if (this.renderable.values.dOrthographic.ref.value !== orthographic) {
      needsUpdate = true;
      ValueCell.update(this.renderable.values.dOrthographic, orthographic);
    }
    ValueCell.updateIfChanged(this.renderable.values.uOutlineThreshold, outlineThreshold);

    if (this.renderable.values.tDepthTransparent.ref.value !== depthTextureTransparent) {
      needsUpdate = true;
      ValueCell.update(this.renderable.values.tDepthTransparent, depthTextureTransparent);
    }
    if (this.renderable.values.tDepthOpaque.ref.value !== depthTextureOpaque) {
      needsUpdate = true;
      ValueCell.update(this.renderable.values.tDepthOpaque, depthTextureOpaque);
    }

    if (needsUpdate) {
      this.renderable.update();
    }

    return { transparentOutline, outlineScale };
  }

  render(b: PassBuilder) {
    b.begin(this.target, { label: 'outlines' });
    b.quad(this.renderable, RenderState.Default);
  }
}

export const OutlinesSchema = {
  ...QuadSchema,
  tDepthOpaque: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  tDepthTransparent: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),

  dOrthographic: DefineSpec('number'),
  uNear: UniformSpec('f'),
  uFar: UniformSpec('f'),
  uInvProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },

  uOutlineThreshold: UniformSpec('f'),
  dTransparentOutline: DefineSpec('boolean'),
};
export const OutlinesShaderCode = ShaderCode('outlines');
export type OutlinesRenderable = ComputeRenderable<Values<typeof OutlinesSchema>>;

export function getOutlinesRenderable(
  ctx: GraphicsContext,
  depthTextureOpaque: Texture,
  depthTextureTransparent: Texture,
  transparentOutline: boolean,
): OutlinesRenderable {
  const width = depthTextureOpaque.getWidth();
  const height = depthTextureOpaque.getHeight();

  const values: Values<typeof OutlinesSchema> = {
    ...QuadValues,
    tDepthOpaque: ValueCell.create(depthTextureOpaque),
    tDepthTransparent: ValueCell.create(depthTextureTransparent),
    uTexSize: ValueCell.create(Vec2.create(width, height)),

    dOrthographic: ValueCell.create(0),
    uNear: ValueCell.create(1),
    uFar: ValueCell.create(10000),
    uInvProjection: ValueCell.create(Mat4.identity()),

    uOutlineThreshold: ValueCell.create(0.33),
    dTransparentOutline: ValueCell.create(transparentOutline),
  };

  const schema = { ...OutlinesSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', OutlinesShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
