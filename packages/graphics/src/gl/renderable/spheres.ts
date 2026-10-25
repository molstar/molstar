/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  type Values,
  InternalSchema,
  SizeSchema,
  type InternalValues,
  ValueSpec,
  DefineSpec,
  GlobalTextureSchema,
  UniformSpec,
  DataSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  InteriorSchema,
  AnimationSchema,
} from './schema.js';
import { SpheresShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const SpheresSchema = {
  ...BaseSchema,
  ...SizeSchema,

  tPositionGroup: DataSpec('rgba32f'),

  padding: ValueSpec('number'),
  uDoubleSided: UniformSpec('b', 'material'),
  dIgnoreLight: DefineSpec('boolean'),
  dCelShaded: DefineSpec('boolean'),
  dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
  dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
  dSolidInterior: DefineSpec('boolean'),
  dClipPrimitive: DefineSpec('boolean'),
  dApproximate: DefineSpec('boolean'),
  uAlphaThickness: UniformSpec('f'),
  uBumpFrequency: UniformSpec('f', 'material'),
  uBumpAmplitude: UniformSpec('f', 'material'),

  lodLevels: ValueSpec('unknown'),
  centerBuffer: ValueSpec('float32'),
  groupBuffer: ValueSpec('float32'),

  ...InteriorSchema,
  ...AnimationSchema,
};
export type SpheresSchema = typeof SpheresSchema;
export type SpheresValues = Values<SpheresSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const SpheresRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...SpheresSchema,
};

export function SpheresRenderable(
  ctx: GraphicsContext,
  id: number,
  values: SpheresValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<SpheresValues> {
  const schema = SpheresRenderableSchema;
  const renderValues: SpheresValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = SpheresShaderCode;
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
