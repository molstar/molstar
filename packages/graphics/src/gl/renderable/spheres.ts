/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  GlobalUniformSchema,
  BaseSchema,
  type Values,
  InternalSchema,
  SizeSchema,
  type InternalValues,
  ValueSpec,
  DefineSpec,
  GlobalTextureSchema,
  UniformSpec,
  DataSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
  InteriorSchema,
  AnimationSchema,
} from './schema.js';
import { SpheresShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra/3d/mat4';
import { type StorageBuffer, StorageBufferSpec } from '../compute.js';

/**
 * The sphere raster's values (WebGPU, `compute/sphere-raster.ts`), set by the renderer per view: in the opaque color
 * pass the impostors leave out the spheres the raster draws (`uRasterMaxRadius` above 0, `rasterizesSphere`); the
 * resolve (the color program with `dSphereResolve`, `RenderItem.renderResolve`) shades the pixels of the object's
 * chunks from the visibility buffer
 */
export const SphereResolveSchema = {
  /** the raster's largest half size in pixels, 0 where it draws none of the object's spheres */
  uRasterMaxRadius: UniformSpec('f'),
  dSphereResolve: DefineSpec('boolean'),
  /** the object's region of the frame's chunks (`SphereRasterTable`) */
  uRasterChunkBase: UniformSpec('i'),
  uRasterChunkCount: UniformSpec('i'),
  /** the visibility buffer's width in pixels */
  uRasterWidth: UniformSpec('i'),
  /** per level (column) its range, fade and scale (`SphereRasterParams.uLevels`) */
  uRasterLevels: UniformSpec('m4'),
  sRasterVisibility: StorageBufferSpec('vec2<u32>', 'read', ['dSphereResolve']),
  sRasterChunks: StorageBufferSpec('vec4<u32>', 'read', ['dSphereResolve']),
  sRasterTransforms: StorageBufferSpec('vec4<f32>', 'read', ['dSphereResolve']),
};
export type SphereResolveValues = Values<typeof SphereResolveSchema>;

/** Values of an object the raster draws none of */
export function createSphereResolveValues(): SphereResolveValues {
  return {
    uRasterMaxRadius: ValueCell.create(0),
    dSphereResolve: ValueCell.create(false),
    uRasterChunkBase: ValueCell.create(0),
    uRasterChunkCount: ValueCell.create(0),
    uRasterWidth: ValueCell.create(0),
    uRasterLevels: ValueCell.create(Mat4()),
    sRasterVisibility: ValueCell.create(undefined as unknown as StorageBuffer),
    sRasterChunks: ValueCell.create(undefined as unknown as StorageBuffer),
    sRasterTransforms: ValueCell.create(undefined as unknown as StorageBuffer),
  };
}

export const SpheresSchema = {
  ...BaseSchema,
  ...SizeSchema,

  tPositionGroup: DataSpec('rgba32f'),

  padding: ValueSpec('number'),
  uDoubleSided: UniformSpec('b', 'material'),
  dIgnoreLight: DefineSpec('boolean'),
  dCelShaded: DefineSpec('boolean'),
  dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
  dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
  dSolidInterior: DefineSpec('boolean'),
  dClipPrimitive: DefineSpec('boolean'),
  dApproximate: DefineSpec('boolean'),
  uAlphaThickness: UniformSpec('f'),
  /** whether the bump noise is compiled (`BaseGeometry.hasBump`; with a substance layer, `dSubstance`, too) */
  dBump: DefineSpec('boolean'),
  uBumpFrequency: UniformSpec('f', 'material'),
  uBumpAmplitude: UniformSpec('f', 'material'),

  lodLevels: ValueSpec('unknown'),

  ...InteriorSchema,
  ...AnimationSchema,
  ...SphereResolveSchema,
};
export type SpheresSchema = typeof SpheresSchema;
export type SpheresValues = Values<SpheresSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const SpheresRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...SpheresSchema,
};

export function SpheresRenderable(
  ctx: GraphicsContext,
  id: number,
  values: SpheresValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<SpheresValues> {
  const schema = SpheresRenderableSchema;
  const renderValues: SpheresValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = SpheresShaderCode;
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
