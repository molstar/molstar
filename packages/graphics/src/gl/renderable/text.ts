/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  UniformSpec,
  type Values,
  InternalSchema,
  SizeSchema,
  type InternalValues,
  TextureSpec,
  DataSpec,
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
  /**
   * Per item 4 texels (`Text.textData`): the label's position and group, depth and whether a tether, corners, texture
   * coordinates; drawn as quads of 4 vertex ids by the shared line indices (`quad-indices.ts`)
   */
  tText: DataSpec('rgba32f'),
  // one font atlas per font, shared by the texts (`.v6/plans/gl/shared-textures.md`)
  tFont: TextureSpec('image-uint8', 'r', 'ubyte', 'linear', undefined, { shared: true }),
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

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const TextRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...TextSchema,
};

export function TextRenderable(
  ctx: GraphicsContext,
  id: number,
  values: TextValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<TextValues> {
  const schema = TextRenderableSchema;
  const renderValues: TextValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = TextShaderCode;
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
