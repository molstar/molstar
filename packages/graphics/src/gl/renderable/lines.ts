/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { WebGLContext } from '../webgl/context.js';
import { createGraphicsRenderItem, type Transparency } from '../webgl/render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  AttributeSpec,
  DefineSpec,
  type Values,
  InternalSchema,
  SizeSchema,
  ElementsSpec,
  type InternalValues,
  GlobalTextureSchema,
  UniformSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  ValueSpec,
  AnimationSchema,
} from './schema.js';
import { ValueCell } from '@molstar/core/util';
import { LinesShaderCode } from '../shader-code.js';

export const LinesSchema = {
  ...BaseSchema,
  ...SizeSchema,
  aGroup: AttributeSpec('float32', 1, 0),
  aMapping: AttributeSpec('float32', 2, 0),
  aStart: AttributeSpec('float32', 3, 0),
  aEnd: AttributeSpec('float32', 3, 0),
  elements: ElementsSpec('uint32'),
  dLineSizeAttenuation: DefineSpec('boolean'),
  uDoubleSided: UniformSpec('b', 'material'),
  dFlipSided: DefineSpec('boolean'),
  stripCount: ValueSpec('number'),
  stripOffsets: ValueSpec('uint32'),

  ...AnimationSchema,
};
export type LinesSchema = typeof LinesSchema;
export type LinesValues = Values<LinesSchema>;

export function LinesRenderable(
  ctx: WebGLContext,
  id: number,
  values: LinesValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<LinesValues> {
  const schema = {
    ...GlobalUniformSchema,
    ...GlobalTextureSchema,
    ...GlobalDefineSchema,
    ...InternalSchema,
    ...LinesSchema,
  };
  const renderValues: LinesValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = LinesShaderCode;
  const renderItem = createGraphicsRenderItem(
    ctx,
    'triangles',
    shaderCode,
    schema,
    renderValues,
    materialId,
    transparency,
  );

  return createRenderable(renderItem, renderValues, state);
}
