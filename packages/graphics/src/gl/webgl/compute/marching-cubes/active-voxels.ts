/**
 * Copyright (c) 2019-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '../../../context.js';
import { type ComputeRenderable, createComputeRenderable } from '../../../renderable.js';
import { RenderState } from '../../../render-state.js';
import { type Values, TextureSpec, UniformSpec, DefineSpec } from '../../../renderable/schema.js';
import type { Texture } from '../../../texture.js';
import { ShaderCode } from '../../../shader-code.js';
import { ValueCell } from '@molstar/core/util';
import type { Vec3, Vec2 } from '@molstar/core/math/linear-algebra';
import { QuadSchema, QuadValues } from '../../../compute/util.js';
import { getTriCount } from './tables.js';
import { isTimingMode } from '@molstar/core/util/debug';
import { setTexelScissor, setTexelViewport } from '../../render-pass-target.js';

/** Channel of the volume data texture holding the value, e.g. `alpha` for densities with a packed group in `rgb` */
export type ValueChannel = 'red' | 'alpha';

export const ActiveVoxelsSchema = {
  ...QuadSchema,

  tTriCount: TextureSpec('image-uint8', 'r', 'ubyte', 'nearest'),
  tVolumeData: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  dValueChannel: DefineSpec('string', ['red', 'alpha']),
  uIsoValue: UniformSpec('f'),

  uGridDim: UniformSpec('v3'),
  uGridTexDim: UniformSpec('v3'),

  uScale: UniformSpec('v2'),
};
export const ActiveVoxelsShaderCode = ShaderCode('active-voxels');
type ActiveVoxelsValues = Values<typeof ActiveVoxelsSchema>;

const ActiveVoxelsName = 'active-voxels';

function getActiveVoxelsRenderable(
  ctx: GraphicsContext,
  volumeData: Texture,
  valueChannel: ValueChannel,
  gridDim: Vec3,
  gridTexDim: Vec3,
  isoValue: number,
  scale: Vec2,
): ComputeRenderable<ActiveVoxelsValues> {
  if (ctx.namedComputeRenderables[ActiveVoxelsName]) {
    const v = ctx.namedComputeRenderables[ActiveVoxelsName].values as ActiveVoxelsValues;

    ValueCell.update(v.uQuadScale, scale);
    ValueCell.update(v.tVolumeData, volumeData);
    ValueCell.update(v.dValueChannel, valueChannel);
    ValueCell.updateIfChanged(v.uIsoValue, isoValue);
    ValueCell.update(v.uGridDim, gridDim);
    ValueCell.update(v.uGridTexDim, gridTexDim);
    ValueCell.update(v.uScale, scale);

    ctx.namedComputeRenderables[ActiveVoxelsName].update();
  } else {
    ctx.namedComputeRenderables[ActiveVoxelsName] = createActiveVoxelsRenderable(
      ctx,
      volumeData,
      valueChannel,
      gridDim,
      gridTexDim,
      isoValue,
      scale,
    );
  }
  return ctx.namedComputeRenderables[ActiveVoxelsName];
}

function createActiveVoxelsRenderable(
  ctx: GraphicsContext,
  volumeData: Texture,
  valueChannel: ValueChannel,
  gridDim: Vec3,
  gridTexDim: Vec3,
  isoValue: number,
  scale: Vec2,
) {
  const values: ActiveVoxelsValues = {
    ...QuadValues,
    tTriCount: ValueCell.create(getTriCount()),

    uQuadScale: ValueCell.create(scale),
    tVolumeData: ValueCell.create(volumeData),
    dValueChannel: ValueCell.create(valueChannel),
    uIsoValue: ValueCell.create(isoValue),
    uGridDim: ValueCell.create(gridDim),
    uGridTexDim: ValueCell.create(gridTexDim),
    uScale: ValueCell.create(scale),
  };

  const schema = { ...ActiveVoxelsSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', ActiveVoxelsShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

export function calcActiveVoxels(
  ctx: GraphicsContext,
  volumeData: Texture,
  valueChannel: ValueChannel,
  gridDim: Vec3,
  gridTexDim: Vec3,
  isoValue: number,
  gridScale: Vec2,
) {
  if (isTimingMode) ctx.timer.group('calcActiveVoxels');
  const { resources } = ctx;
  const width = volumeData.getWidth();
  const height = volumeData.getHeight();

  if (!ctx.namedFramebuffers[ActiveVoxelsName]) {
    ctx.namedFramebuffers[ActiveVoxelsName] = resources.framebuffer();
  }
  const framebuffer = ctx.namedFramebuffers[ActiveVoxelsName];

  if (!ctx.namedTextures[ActiveVoxelsName]) {
    ctx.namedTextures[ActiveVoxelsName] = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
  }
  const activeVoxelsTex = ctx.namedTextures[ActiveVoxelsName];
  activeVoxelsTex.define(width, height);

  const renderable = getActiveVoxelsRenderable(ctx, volumeData, valueChannel, gridDim, gridTexDim, isoValue, gridScale);

  activeVoxelsTex.attachFramebuffer(framebuffer, 0);
  ctx.beginRenderPass(framebuffer, { label: 'active voxels' });
  ctx.setRenderState(RenderState.Default);
  setTexelViewport(ctx, 0, 0, width, height);
  setTexelScissor(ctx, 0, 0, width, height);
  ctx.clear({ color: [0, 0, 0, 0] });
  setTexelScissor(ctx, 0, 0, gridTexDim[0], gridTexDim[1]);
  renderable.render();

  // console.log('gridScale', gridScale, 'gridTexDim', gridTexDim, 'gridDim', gridDim);
  // console.log('volumeData', volumeData);
  // console.log('at', readTexture(ctx, activeVoxelsTex));
  // printTextureImage(readTexture(ctx, activeVoxelsTex), { scale: 0.75 });

  if (isTimingMode) ctx.timer.groupEnd('calcActiveVoxels');

  return activeVoxelsTex;
}
