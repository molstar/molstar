/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type ComputeRenderable, createComputeRenderable } from '../../../renderable.js';
import type { WebGLContext } from '../../context.js';
import { RenderState } from '../../../render-state.js';
import { type Values, TextureSpec } from '../../../renderable/schema.js';
import type { Texture } from '../../../texture.js';
import { ShaderCode } from '../../../shader-code.js';
import { ValueCell } from '@molstar/core/util';
import { unpackRGBToInt } from '@molstar/core/util/number-packing';
import { QuadSchema, QuadValues } from '../../../compute/util.js';
import { isTimingMode } from '@molstar/core/util/debug';
import { setTexelViewport } from '../../render-pass-target.js';
import { createScratch } from '../../../scratch.js';
import type { Framebuffer } from '../../framebuffer.js';
import { asGLTexture } from '../../texture.js';

/** the scratch resources of this module, per context (`gl/scratch.ts`) */
const ScratchTextures = createScratch<Texture>('histogram-pyramid/sum textures');
const ScratchFramebuffers = createScratch<Framebuffer>('histogram-pyramid/sum framebuffers');
const ScratchRenderables = createScratch<ComputeRenderable<any>>('histogram-pyramid/sum renderables');

export const HistopyramidSumSchema = {
  ...QuadSchema,
  tTexture: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, { sampleType: 'sint' }),
};
export const HistopyramidSumShaderCode = ShaderCode('sum');
type HistopyramidSumValues = Values<typeof HistopyramidSumSchema>;

const HistopyramidSumName = 'histopyramid-sum';

function getHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture): ComputeRenderable<HistopyramidSumValues> {
  if (ScratchRenderables.get(ctx, HistopyramidSumName)) {
    const v = ScratchRenderables.get(ctx, HistopyramidSumName)!.values as HistopyramidSumValues;

    ValueCell.update(v.tTexture, texture);

    ScratchRenderables.get(ctx, HistopyramidSumName)!.update();
  } else {
    ScratchRenderables.set(ctx, createHistopyramidSumRenderable(ctx, texture), HistopyramidSumName);
  }
  return ScratchRenderables.get(ctx, HistopyramidSumName)!;
}

function createHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture) {
  const values: HistopyramidSumValues = {
    ...QuadValues,
    tTexture: ValueCell.create(texture),
  };

  const schema = { ...HistopyramidSumSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', HistopyramidSumShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

const sumBytes = new Uint8Array(4);
const sumInts = new Int32Array(4);

export function getHistopyramidSum(ctx: WebGLContext, pyramidTopTexture: Texture) {
  if (isTimingMode) ctx.timer.group('getHistopyramidSum');
  const { gl, resources, capabilities } = ctx;

  const renderable = getHistopyramidSumRenderable(ctx, pyramidTopTexture);

  const framebuffer = ScratchFramebuffers.getOrCreate(ctx, () => resources.framebuffer(), HistopyramidSumName);

  if (!ScratchTextures.get(ctx, HistopyramidSumName)) {
    ScratchTextures.set(
      ctx,
      capabilities.integerTextures
        ? resources.texture('image-int32', 'rgba', 'int', 'nearest')
        : resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest'),
      HistopyramidSumName,
    );
    ScratchTextures.get(ctx, HistopyramidSumName)!.define(1, 1);
  }
  const sumTexture = ScratchTextures.get(ctx, HistopyramidSumName)!;
  asGLTexture(sumTexture).attachFramebuffer(framebuffer, 0);

  ctx.beginRenderPass(framebuffer, { label: 'histogram pyramid sum' });
  ctx.setRenderState(RenderState.Default);
  ctx.disableScissor();
  setTexelViewport(ctx, 0, 0, 1, 1);
  renderable.render();

  ctx.readPixels(0, 0, 1, 1, capabilities.integerTextures ? sumInts : sumBytes);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  if (isTimingMode) ctx.timer.groupEnd('getHistopyramidSum');

  return capabilities.integerTextures ? sumInts[0] : unpackRGBToInt(sumBytes[0], sumBytes[1], sumBytes[2]);
}
