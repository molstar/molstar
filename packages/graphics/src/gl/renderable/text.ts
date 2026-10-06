/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
  UniformSpec,
  type Values,
  InternalSchema,
  SizeSchema,
  type InternalValues,
  TextureSpec,
  ElementsSpec,
  ValueSpec,
  GlobalTextureSchema,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
} from './schema.js';
import { TextShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const TextSchema = {
  ...BaseSchema,
  ...SizeSchema,
  aGroup: AttributeSpec('float32', 1, 0),
  aPosition: AttributeSpec('float32', 3, 0),
  aMapping: AttributeSpec('float32', 2, 0),
  aDepth: AttributeSpec('float32', 1, 0),
  elements: ElementsSpec('uint32'),

  aTexCoord: AttributeSpec('float32', 2, 0),
  tFont: TextureSpec('image-uint8', 'alpha', 'ubyte', 'linear'),
  padding: ValueSpec('number'),

  uBorderWidth: UniformSpec('f', 'material'),
  uBorderColor: UniformSpec('v3', 'material'),
  uOffsetX: UniformSpec('f', 'material'),
  uOffsetY: UniformSpec('f', 'material'),
  uOffsetZ: UniformSpec('f', 'material'),
  uBackgroundColor: UniformSpec('v3', 'material'),
  uBackgroundOpacity: UniformSpec('f', 'material'),
};
export type TextSchema = typeof TextSchema;
export type TextValues = Values<TextSchema>;

export function TextRenderable(
  ctx: WebGLContext,
  id: number,
  values: TextValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<TextValues> {
  const schema = {
    ...GlobalUniformSchema,
    ...GlobalTextureSchema,
    ...GlobalDefineSchema,
    ...InternalSchema,
    ...TextSchema,
  };
  const renderValues: TextValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = TextShaderCode;
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
