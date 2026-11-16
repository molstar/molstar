/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createInstanceStorage, usesInstanceStorage, withInstanceStorage } from '../instance-storage.js';
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
import { Vec2 } from '@molstar/core/math/linear-algebra/3d/vec2';
import { type StorageBuffer, StorageBufferSpec } from '../compute.js';
import { createDataArray } from './util.js';

/**
 * The spheres' second color, blended in near the camera (a particle target's sites by element):
 * by a sphere's radius on screen, from the primary color at `uSiteColorRange.x` pixels to
 * the second at `uSiteColorRange.y`
 */
export const SiteColorSchema = {
  /**
   * off; by the site's type (`dSites`: a texel a type, a merged object's the members' tables packed, a segment's first
   * entry in its data); by group (a texel a site)
   */
  dSiteColor: DefineSpec('string', ['off', 'siteType', 'group']),
  /** rgb the second color, alpha 0 where it is the primary color (a carbon the particle's) */
  tSiteColor: DataSpec('rgba8', ['dSiteColor_siteType', 'dSiteColor_group']),
  /** the radius on screen in pixels from which the second color blends in (x) and from which it is whole (y) */
  uSiteColorRange: UniformSpec('v2', 'material'),
};
export type SiteColorValues = Values<typeof SiteColorSchema>;

/** Values of spheres without a second color */
export function createSiteColorValues(): SiteColorValues {
  return {
    dSiteColor: ValueCell.create('off'),
    tSiteColor: ValueCell.create(createDataArray(1, 4, Uint8Array)),
    uSiteColorRange: ValueCell.create(Vec2()),
  };
}

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
  /** the view's origin in its target in pixels: the visibility buffer's pixel of a fragment its coordinates less it */
  uRasterOrigin: UniformSpec('v2'),
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
    uRasterOrigin: ValueCell.create(Vec2()),
    uRasterLevels: ValueCell.create(Mat4()),
    sRasterVisibility: ValueCell.create(undefined as unknown as StorageBuffer),
    sRasterChunks: ValueCell.create(undefined as unknown as StorageBuffer),
    sRasterTransforms: ValueCell.create(undefined as unknown as StorageBuffer),
  };
}

export const SpheresSchema = {
  ...BaseSchema,
  ...SizeSchema,

  tPositionGroup: DataSpec('rgb32f-i32', ['!dSites']),
  /** quantized centers (`dSites`, `SphereSites`): int16 x, y, z a sphere, its center `uSiteOrigin + q * uSiteScale`; group the sphere */
  tSite: DataSpec('rgba16i', ['dSites']),
  uSiteOrigin: UniformSpec('v3'),
  uSiteScale: UniformSpec('f'),
  dSites: DefineSpec('boolean'),

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
  /** the stride of the level drawn (per level list, `lodLevels`): quad `p` draws item `p * uLodStride`; 1 without */
  uLodStride: UniformSpec('f', undefined, 'draw'),

  ...InteriorSchema,
  ...AnimationSchema,
  ...SphereResolveSchema,
  ...SiteColorSchema,
};
export type SpheresSchema = typeof SpheresSchema;
export type SpheresValues = Values<SpheresSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const SpheresRenderableSchema = withInstanceStorage({
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...SpheresSchema,
});

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
  // poses read by index from storage where the context does (`instance-storage.ts`), not packed as rows
  const storage = usesInstanceStorage(ctx, values) ? createInstanceStorage(ctx.compute!, renderValues) : undefined;
  Object.assign(renderValues, storage?.values ?? { dInstanceStorage: ValueCell.create(false) });
  const shaderCode = SpheresShaderCode;
  const renderItem = ctx.createGraphicsRenderItem(
    'triangles',
    shaderCode,
    schema,
    renderValues,
    materialId,
    transparency,
  );
  return createRenderable(renderItem, renderValues, state, storage);
}
