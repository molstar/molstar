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
  type Values,
  InternalSchema,
  SizeSchema,
  DataSpec,
  type InternalValues,
  GlobalTextureSchema,
  UniformSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  AnimationSchema,
} from './schema.js';
import { ValueCell } from '@molstar/core/util';
import { LinesShaderCode } from '../shader-code.js';

export const LinesSchema = {
  ...BaseSchema,
  ...SizeSchema,
  /**
   * Per line 2 texels (`Lines.lineData`): start xyz and group, end xyz and whether it continues the line before (strips);
   * drawn as quads of 4 vertex ids by the shared line indices (`quad-indices.ts`)
   */
  tLine: DataSpec('rgba32f'),
  dLineSizeAttenuation: DefineSpec('boolean'),
  uDoubleSided: UniformSpec('b', 'material'),
  dFlipSided: DefineSpec('boolean'),

  ...AnimationSchema,
};
export type LinesSchema = typeof LinesSchema;
export type LinesValues = Values<LinesSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const LinesRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...LinesSchema,
};

export function LinesRenderable(
  ctx: GraphicsContext,
  id: number,
  values: LinesValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<LinesValues> {
  const schema = LinesRenderableSchema;
  const renderValues: LinesValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = LinesShaderCode;
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
