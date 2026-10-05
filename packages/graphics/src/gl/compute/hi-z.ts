/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { WebGLContext } from '../webgl/context.js';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { ValueCell } from '@molstar/core/util/value-cell';
import { type ComputeRenderable, createComputeRenderable } from '../renderable.js';
import { TextureSpec, UniformSpec, type Values } from '../renderable/schema.js';
import { ShaderCode } from '../shader-code.js';
import { hiZ_frag } from '../shader/hi-z.frag.js';
import { quad_vert } from '../shader/quad.vert.js';
import { createComputeRenderItem } from '../webgl/render-item.js';
import type { Texture } from '../webgl/texture.js';
import { QuadSchema, QuadValues } from './util.js';

const HiZSchema = {
  ...QuadSchema,
  tPreviousLevel: TextureSpec('texture', 'alpha', 'float', 'nearest'),
  uInvSize: UniformSpec('v2'),
  uOffset: UniformSpec('v2'),
};
const HiZShaderCode = ShaderCode('hi-z', quad_vert, hiZ_frag);
export type HiZRenderable = ComputeRenderable<Values<typeof HiZSchema>>;

export function createHiZRenderable(ctx: WebGLContext, previousLevel: Texture): HiZRenderable {
  const values: Values<typeof HiZSchema> = {
    ...QuadValues,
    tPreviousLevel: ValueCell.create(previousLevel),
    uInvSize: ValueCell.create(Vec2()),
    uOffset: ValueCell.create(Vec2()),
  };

  const schema = { ...HiZSchema };
  const renderItem = createComputeRenderItem(ctx, 'triangles', HiZShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}
