/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
  AttributeSpec,
  ElementsSpec,
  DefineSpec,
  type Values,
  InternalSchema,
  type InternalValues,
  GlobalTextureSchema,
  ValueSpec,
  UniformSpec,
  GlobalDefineSchema,
  type GlobalDefineValues,
  type GlobalDefines,
  InteriorSchema,
  AnimationSchema,
} from './schema.js';
import { MeshShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const MeshSchema = {
  ...BaseSchema,
  aGroup: AttributeSpec('float32', 1, 0),
  aPosition: AttributeSpec('float32', 3, 0),
  /** unit normals, three floats on the CPU, octahedral in two snorm16 on the GPU (`oct-normals.ts`) */
  aNormal: AttributeSpec('oct16', 2, 0),
  elements: ElementsSpec('uint32'),
  dVaryingGroup: DefineSpec('boolean'),
  dFlatShaded: DefineSpec('boolean'),
  uDoubleSided: UniformSpec('b', 'material'),
  dFlipSided: DefineSpec('boolean'),
  dIgnoreLight: DefineSpec('boolean'),
  dCelShaded: DefineSpec('boolean'),
  dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
  dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
  dSolidInterior: DefineSpec('boolean'),
  /** whether the bump noise is compiled (`BaseGeometry.hasBump`; with a substance layer, `dSubstance`, too) */
  dBump: DefineSpec('boolean'),
  uBumpFrequency: UniformSpec('f', 'material'),
  uBumpAmplitude: UniformSpec('f', 'material'),
  /** Å toward the camera, in depth only (`Mesh.Params.depthOffset`) */
  uDepthOffset: UniformSpec('f'),
  /**
   * CSS pixels per unit of the mesh, constant on screen: the vertex stage scales it about the instance's origin by the
   * world size of a pixel there (the instance's own scale multiplies); 0: world units
   */
  uScreenSize: UniformSpec('f'),
  meta: ValueSpec('unknown'),

  ...InteriorSchema,
  ...AnimationSchema,
} as const;
export type MeshSchema = typeof MeshSchema;
export type MeshValues = Values<MeshSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const MeshRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...MeshSchema,
};

export function MeshRenderable(
  ctx: GraphicsContext,
  id: number,
  values: MeshValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<MeshValues> {
  const schema = MeshRenderableSchema;
  const renderValues: MeshValues & InternalValues & GlobalDefineValues = {
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

  return createRenderable(renderItem, renderValues, state);
}
