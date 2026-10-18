/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '../context.js';
import { createNullTexture, type Texture } from '../texture.js';
import { ValueCell } from '@molstar/core/util';
import { ValueSpec, AttributeSpec, UniformSpec, type Values, TextureSpec, DefineSpec } from '../renderable/schema.js';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { ShaderCode } from '../shader-code.js';
import { type ComputeRenderable, createComputeRenderable } from '../renderable.js';

export const QuadPositions = new Float32Array([
  1.0,
  1.0,
  -1.0,
  1.0,
  -1.0,
  -1.0, // First triangle
  -1.0,
  -1.0,
  1.0,
  -1.0,
  1.0,
  1.0, // Second triangle
]);

export const QuadSchema = {
  drawCount: ValueSpec('number'),
  instanceCount: ValueSpec('number'),
  aPosition: AttributeSpec('float32', 2, 0),
  uQuadScale: UniformSpec('v2'),
};

export const QuadValues: Values<typeof QuadSchema> = {
  drawCount: ValueCell.create(6),
  instanceCount: ValueCell.create(1),
  aPosition: ValueCell.create(QuadPositions),
  uQuadScale: ValueCell.create(Vec2.create(1, 1)),
};

//

export const CopySchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),
  /** pack the red channel (depth) into RGBA, for depth copies into non-float targets */
  dPackDepth: DefineSpec('boolean'),
};
export const CopyShaderCode = ShaderCode('copy');
export type CopyRenderable = ComputeRenderable<Values<typeof CopySchema>>;

export function createCopyRenderable(ctx: GraphicsContext, texture: Texture, packDepth = false): CopyRenderable {
  const values: Values<typeof CopySchema> = {
    ...QuadValues,
    tColor: ValueCell.create(texture),
    uTexSize: ValueCell.create(Vec2.create(texture.getWidth(), texture.getHeight())),
    dPackDepth: ValueCell.create(packDepth),
  };

  const schema = { ...CopySchema };
  const renderItem = ctx.createComputeRenderItem('triangles', CopyShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

const SharedCopyName = 'shared-copy';

export function getSharedCopyRenderable(ctx: GraphicsContext, texture: Texture) {
  if (!ctx.namedComputeRenderables[SharedCopyName]) {
    ctx.namedComputeRenderables[SharedCopyName] = createCopyRenderable(ctx, createNullTexture());
  }
  const copy = ctx.namedComputeRenderables[SharedCopyName] as CopyRenderable;
  ValueCell.update(copy.values.tColor, texture);
  ValueCell.update(
    copy.values.uTexSize,
    Vec2.set(copy.values.uTexSize.ref.value, texture.getWidth(), texture.getHeight()),
  );
  copy.update();
  return copy;
}
