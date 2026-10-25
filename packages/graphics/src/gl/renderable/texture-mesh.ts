/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
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
  type InternalValues,
  UniformSpec,
  TextureSpec,
  GlobalTextureSchema,
  ValueSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  InteriorSchema,
  AnimationSchema,
} from './schema.js';
import { MeshShaderCode } from '../shader-code.js';
import { StorageBufferSpec } from '../compute.js';
import { ValueCell } from '@molstar/core/util';

export const TextureMeshSchema = {
  ...BaseSchema,
  uGeoTexDim: UniformSpec('v2', 'buffered'),
  tPosition: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tGroup: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tNormal: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  /** vertices in pages of the vertex arena instead of the textures (WebGPU), see `computeChunkedIsosurface` */
  dIsoChunks: DefineSpec('boolean'),
  sIsoTable: StorageBufferSpec('u32', 'read', ['dIsoChunks']),
  sIsoVertices: StorageBufferSpec('u32', 'read', ['dIsoChunks']),
  dVaryingGroup: DefineSpec('boolean'),
  dFlatShaded: DefineSpec('boolean'),
  uDoubleSided: UniformSpec('b', 'material'),
  dFlipSided: DefineSpec('boolean'),
  dIgnoreLight: DefineSpec('boolean'),
  dCelShaded: DefineSpec('boolean'),
  dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
  dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
  dSolidInterior: DefineSpec('boolean'),
  uBumpFrequency: UniformSpec('f', 'material'),
  uBumpAmplitude: UniformSpec('f', 'material'),
  meta: ValueSpec('unknown'),

  ...InteriorSchema,
  ...AnimationSchema,
};
export type TextureMeshSchema = typeof TextureMeshSchema;
export type TextureMeshValues = Values<TextureMeshSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const TextureMeshRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...TextureMeshSchema,
};

export function TextureMeshRenderable(
  ctx: GraphicsContext,
  id: number,
  values: TextureMeshValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<TextureMeshValues> {
  const schema = TextureMeshRenderableSchema;
  const renderValues: TextureMeshValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = MeshShaderCode;
  const renderItem = ctx.createGraphicsRenderItem(
    'triangles',
    shaderCode,
    schema,
    renderValues,
    materialId,
    transparency,
  );

  const renderable = createRenderable(renderItem, renderValues, state);
  const getItemByteCount = renderable.getByteCount;
  // with its pages of the shared vertex arena (see `computeChunkedIsosurface`), known once read back
  return Object.assign(renderable, {
    getByteCount: () => {
      const meta = values.meta.ref.value as { chunkByteCount?: number } | undefined;
      return getItemByteCount() + (values.dIsoChunks.ref.value ? (meta?.chunkByteCount ?? 0) : 0);
    },
  });
}
