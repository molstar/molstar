/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  DefineSpec,
  DataSpec,
  ValueSpec,
  type Values,
  InternalSchema,
  SizeSchema,
  type InternalValues,
  GlobalTextureSchema,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  AnimationSchema,
} from './schema.js';
import { PointsShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const PointsSchema = {
  ...BaseSchema,
  ...SizeSchema,
  tPositionGroup: DataSpec('rgba32f'),
  centerBuffer: ValueSpec('float32'),
  groupBuffer: ValueSpec('float32'),
  dPointSizeAttenuation: DefineSpec('boolean'),
  dPointStyle: DefineSpec('string', ['square', 'circle', 'fuzzy']),

  ...AnimationSchema,
};
export type PointsSchema = typeof PointsSchema;
export type PointsValues = Values<PointsSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const PointsRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...PointsSchema,
};

export function PointsRenderable(
  ctx: GraphicsContext,
  id: number,
  values: PointsValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<PointsValues> {
  const schema = PointsRenderableSchema;
  const renderValues: PointsValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = PointsShaderCode;
  // each point is a quad of two triangles, see the points vertex shader
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
