/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '../context.js';
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { Vec4 } from '@molstar/core/math/linear-algebra/3d/vec4';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { ValueCell } from '@molstar/core/util/value-cell';
import { Clip } from '@molstar/core/util/clip';
import { ClipStorageVec4s } from '../clip-storage.js';
import { type ComputeRenderable, createComputeRenderable } from '../renderable.js';
import { TextureSpec, UniformSpec, type Values } from '../renderable/schema.js';
import { ShaderCode } from '../shader-code.js';
import type { Texture } from '../texture.js';
import { QuadSchema, QuadValues } from './util.js';

export const HiZSchema = {
  ...QuadSchema,
  tPreviousLevel: TextureSpec('texture', 'r', 'float', 'nearest'),
  uInvSize: UniformSpec('v2'),
  /** the region of the previous level sampled, in uv (`xy` its corner, `zw` its extent) */
  uRect: { ...UniformSpec('v4'), windowRect: 'rect' as const },
  /** size of the level rendered */
  uTexSize: UniformSpec('v2'),
};
export const HiZShaderCode = ShaderCode('hi-z');
export type HiZRenderable = ComputeRenderable<Values<typeof HiZSchema>>;

export function createHiZRenderable(ctx: GraphicsContext, previousLevel: Texture): HiZRenderable {
  const values: Values<typeof HiZSchema> = {
    ...QuadValues,
    tPreviousLevel: ValueCell.create(previousLevel),
    uInvSize: ValueCell.create(Vec2()),
    uRect: ValueCell.create(Vec4()),
    uTexSize: ValueCell.create(Vec2()),
  };

  const schema = { ...HiZSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', HiZShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

/**
 * The first level primed from an occluder grid (WebGPU, `OccluderGrid`): the depth texture's farthest as `HiZSchema`'s,
 * the nearer where a ray through the grid's optical densities reaches `uOpticalDepth` (the gaps between molecules
 * closed); the grid's cells in a 2d texture, its slices along z side by side (`uGridColumns` a row)
 */
export const HiZOccluderSchema = {
  ...HiZSchema,
  tGrid: TextureSpec('texture', 'r', 'float', 'nearest'),
  /** the viewport in the drawing buffer's uv (top-left origin): x, y, width, height */
  uViewportRect: UniformSpec('v4'),
  /** the camera's, mapped by the backend to its reversed clip depth (`DepthRangeTag`): near at 1 */
  uInvViewProjection: { ...UniformSpec('m4'), depthRange: 'inverse' as const },
  uViewProjection: { ...UniformSpec('m4'), depthRange: 'projection' as const },
  uGridMin: UniformSpec('v3'),
  uGridMax: UniformSpec('v3'),
  uGridDims: UniformSpec('v3'),
  uGridColumns: UniformSpec('f'),
  uOpticalDepth: UniformSpec('f'),
  /**
   * the grid's object's clip objects (`clip_storage_wgsl`: `ClipStorageVec4s` vec4s each, `uClipCount` of them): what
   * they remove hides nothing (a ray's samples there skipped), as the object draws it; its uniforms each frame, a cut
   * moved (a drag) without the grid uploaded again
   */
  uClipCount: UniformSpec('i'),
  uClipData: UniformSpec('v4[]', undefined, 'object', Clip.MaxObjectCount * ClipStorageVec4s),
};
export const HiZOccluderShaderCode = ShaderCode('hi-z-occluder');
export type HiZOccluderRenderable = ComputeRenderable<Values<typeof HiZOccluderSchema>>;

export function createHiZOccluderRenderable(
  ctx: GraphicsContext,
  previousLevel: Texture,
  grid: Texture,
): HiZOccluderRenderable {
  const values: Values<typeof HiZOccluderSchema> = {
    ...QuadValues,
    tPreviousLevel: ValueCell.create(previousLevel),
    uInvSize: ValueCell.create(Vec2()),
    uRect: ValueCell.create(Vec4()),
    uTexSize: ValueCell.create(Vec2()),
    tGrid: ValueCell.create(grid),
    uViewportRect: ValueCell.create(Vec4()),
    uInvViewProjection: ValueCell.create(Mat4()),
    uViewProjection: ValueCell.create(Mat4()),
    uGridMin: ValueCell.create(Vec3()),
    uGridMax: ValueCell.create(Vec3()),
    uGridDims: ValueCell.create(Vec3()),
    uGridColumns: ValueCell.create(1),
    uOpticalDepth: ValueCell.create(7),
    uClipCount: ValueCell.create(0),
    uClipData: ValueCell.create(new Array<number>(Clip.MaxObjectCount * ClipStorageVec4s * 4).fill(0)),
  };

  const renderItem = ctx.createComputeRenderItem('triangles', HiZOccluderShaderCode, { ...HiZOccluderSchema }, values);
  return createComputeRenderable(renderItem, values);
}
