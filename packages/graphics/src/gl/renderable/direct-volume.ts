/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { GraphicsContext } from '../context.js';
import type { Transparency } from '../render-item.js';
import {
  AttributeSpec,
  type Values,
  UniformSpec,
  GlobalUniformSchema,
  InternalSchema,
  TextureSpec,
  ElementsSpec,
  DefineSpec,
  type InternalValues,
  GlobalTextureSchema,
  BaseSchema,
  ValueSpec,
  type GlobalDefineValues,
  type GlobalDefines,
  GlobalDefineSchema,
} from './schema.js';
import { DirectVolumeShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const DirectVolumeSchema = {
  ...BaseSchema,

  aPosition: AttributeSpec('float32', 3, 0),
  elements: ElementsSpec('uint32'),

  uBboxMin: UniformSpec('v3'),
  uBboxMax: UniformSpec('v3'),
  uBboxSize: UniformSpec('v3'),
  uMaxSteps: { ...UniformSpec('i'), loopBound: true },
  uStepScale: UniformSpec('f'),
  uJumpLength: UniformSpec('f'),
  uTransform: UniformSpec('m4'),
  uGridDim: UniformSpec('v3'),
  tTransferTex: TextureSpec('image-uint8', 'r', 'ubyte', 'linear'),
  /** the transfer function's opacity between two values (`preIntegration`, `createPreIntegratedTransferTexture`) */
  tPreIntegratedTex: TextureSpec('image-uint8', 'r', 'ubyte', 'linear', undefined, { usedIf: ['dPreIntegration'] }),
  dPreIntegration: DefineSpec('boolean'),
  uTransferScale: UniformSpec('f', 'material'),
  /** the opacity at which the marking and pick variants put the volume's surface */
  uSurfaceOpacity: UniformSpec('f'),
  /** `renderMode`: composite blends the samples of a ray, maximum keeps its brightest one, attenuated its brightest
   * dimmed by the density in front, additive sums them */
  dRenderMode: DefineSpec('string', ['composite', 'maximum', 'attenuated', 'additive']),
  /** the value the `maximum` mode's window starts at: a sample below it adds nothing */
  uMaximumStart: UniformSpec('f'),
  /** `attenuation`: the dimming of the `attenuated` mode, the gain of the `additive` one */
  uAttenuation: UniformSpec('f'),
  /** the grid's bricks (`BrickRange`), their visibility under the transfer function per brick and the bricks per axis */
  brickRange: ValueSpec('unknown'),
  tBrickTex: TextureSpec('image-uint8', 'r', 'ubyte', 'nearest'),
  uBrickTexDim: UniformSpec('v2'),
  uBrickDim: UniformSpec('v3'),

  dGridTexType: DefineSpec('string', ['2d', '3d']),
  uGridTexDim: UniformSpec('v3'),
  tGridTex: TextureSpec('texture', 'rgba', 'ubyte', 'linear', undefined, { dimension: { define: 'dGridTexType' } }),

  uCellDim: UniformSpec('v3'),
  uCartnToUnit: UniformSpec('m4'),
  uUnitToCartn: UniformSpec('m4'),
  dPackedGroup: DefineSpec('boolean'),
  dAxisOrder: DefineSpec('string', ['012', '021', '102', '120', '201', '210']),

  dIgnoreLight: DefineSpec('boolean'),
  dCelShaded: DefineSpec('boolean'),
  dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),

  meta: ValueSpec('unknown'),
};
export type DirectVolumeSchema = typeof DirectVolumeSchema;
export type DirectVolumeValues = Values<DirectVolumeSchema>;

/** The schema of the renderable: the object's with the global and internal ones (one object, shared by the items) */
export const DirectVolumeRenderableSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...DirectVolumeSchema,
};

export function DirectVolumeRenderable(
  ctx: GraphicsContext,
  id: number,
  values: DirectVolumeValues,
  state: RenderableState,
  materialId: number,
  transparency: Transparency,
  globals: GlobalDefines,
): Renderable<DirectVolumeValues> {
  const schema = DirectVolumeRenderableSchema;
  const renderValues: DirectVolumeValues & InternalValues & GlobalDefineValues = {
    ...values,
    uObjectId: ValueCell.create(id),
    dLightCount: ValueCell.create(globals.dLightCount),
    dColorMarker: ValueCell.create(globals.dColorMarker),
  };
  const shaderCode = DirectVolumeShaderCode;
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
