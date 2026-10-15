/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  DefineSpec,
  type Values,
  InternalSchema,
  type InternalValues,
  UniformSpec,
  TextureSpec,
  GlobalTextureSchema,
  ValueSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  InteriorSchema,
  AnimationSchema,
} from './schema.js';
import { MeshShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const TextureMeshSchema = {
  ...BaseSchema,
  uGeoTexDim: UniformSpec('v2', 'buffered'),
  tPosition: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tGroup: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tNormal: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  dVaryingGroup: DefineSpec('boolean'),
  dFlatShaded: DefineSpec('boolean'),
  uDoubleSided: UniformSpec('b', 'material'),
  dFlipSided: DefineSpec('boolean'),
  dIgnoreLight: DefineSpec('boolean'),
  dCelShaded: DefineSpec('boolean'),
  dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
  dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
  dSolidInterior: DefineSpec('boolean'),
  uBumpFrequency: UniformSpec('f', 'material'),
  uBumpAmplitude: UniformSpec('f', 'material'),
  meta: ValueSpec('unknown'),

  ...InteriorSchema,
  ...AnimationSchema,
};
export type TextureMeshSchema = typeof TextureMeshSchema;
export type TextureMeshValues = Values<TextureMeshSchema>;

export function TextureMeshRenderable(
  ctx: GraphicsContext,
  id: number,
  values: TextureMeshValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<TextureMeshValues> {
  const schema = {
    ...GlobalUniformSchema,
    ...GlobalTextureSchema,
    ...GlobalDefineSchema,
    ...InternalSchema,
    ...TextureMeshSchema,
  };
  const renderValues: TextureMeshValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = MeshShaderCode;
  const renderItem = ctx.createGraphicsRenderItem(
    'triangles',
    shaderCode,
    schema,
    renderValues,
    materialId,
    transparency,
  );

  return createRenderable(renderItem, renderValues, state);
}
