/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '../context.js';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { ValueCell } from '@molstar/core/util/value-cell';
import { type ComputeRenderable, createComputeRenderable } from '../renderable.js';
import { TextureSpec, UniformSpec, type Values } from '../renderable/schema.js';
import { ShaderCode } from '../shader-code.js';
import type { Texture } from '../texture.js';
import { QuadSchema, QuadValues } from './util.js';

export const HiZSchema = {
  ...QuadSchema,
  tPreviousLevel: TextureSpec('texture', 'r', 'float', 'nearest'),
  uInvSize: UniformSpec('v2'),
  uOffset: UniformSpec('v2'),
};
export const HiZShaderCode = ShaderCode('hi-z');
export type HiZRenderable = ComputeRenderable<Values<typeof HiZSchema>>;

export function createHiZRenderable(ctx: GraphicsContext, previousLevel: Texture): HiZRenderable {
  const values: Values<typeof HiZSchema> = {
    ...QuadValues,
    tPreviousLevel: ValueCell.create(previousLevel),
    uInvSize: ValueCell.create(Vec2()),
    uOffset: ValueCell.create(Vec2()),
  };

  const schema = { ...HiZSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', HiZShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
