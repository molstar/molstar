/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { WebGLContext } from '../../context.js';
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
import { createScratch } from '../../../scratch.js';
import type { Framebuffer } from '../../framebuffer.js';
import { asGLTexture } from '../../texture.js';

/** the scratch resources of this module, per context (`gl/scratch.ts`) */
const ScratchTextures = createScratch<Texture>('marching-cubes/active-voxels textures');
const ScratchFramebuffers = createScratch<Framebuffer>('marching-cubes/active-voxels framebuffers');
const ScratchRenderables = createScratch<ComputeRenderable<any>>('marching-cubes/active-voxels renderables');

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
  ctx: WebGLContext,
  volumeData: Texture,
  valueChannel: ValueChannel,
  gridDim: Vec3,
  gridTexDim: Vec3,
  isoValue: number,
  scale: Vec2,
): ComputeRenderable<ActiveVoxelsValues> {
  if (ScratchRenderables.get(ctx, ActiveVoxelsName)) {
    const v = ScratchRenderables.get(ctx, ActiveVoxelsName)!.values as ActiveVoxelsValues;

    ValueCell.update(v.uQuadScale, scale);
    ValueCell.update(v.tVolumeData, volumeData);
    ValueCell.update(v.dValueChannel, valueChannel);
    ValueCell.updateIfChanged(v.uIsoValue, isoValue);
    ValueCell.update(v.uGridDim, gridDim);
    ValueCell.update(v.uGridTexDim, gridTexDim);
    ValueCell.update(v.uScale, scale);

    ScratchRenderables.get(ctx, ActiveVoxelsName)!.update();
  } else {
    ScratchRenderables.set(
      ctx,
      createActiveVoxelsRenderable(ctx, volumeData, valueChannel, gridDim, gridTexDim, isoValue, scale),
      ActiveVoxelsName,
    );
  }
  return ScratchRenderables.get(ctx, ActiveVoxelsName)!;
}

function createActiveVoxelsRenderable(
  ctx: WebGLContext,
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
  ctx: WebGLContext,
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

  const framebuffer = ScratchFramebuffers.getOrCreate(ctx, () => resources.framebuffer(), ActiveVoxelsName);

  const activeVoxelsTex = ScratchTextures.getOrCreate(
    ctx,
    () => resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest'),
    ActiveVoxelsName,
  );
  activeVoxelsTex.define(width, height);

  const renderable = getActiveVoxelsRenderable(ctx, volumeData, valueChannel, gridDim, gridTexDim, isoValue, gridScale);

  asGLTexture(activeVoxelsTex).attachFramebuffer(framebuffer, 0);
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
