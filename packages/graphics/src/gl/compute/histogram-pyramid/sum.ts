/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type ComputeRenderable, createComputeRenderable } from '../../renderable.js';
import type { WebGLContext } from '../../webgl/context.js';
import { RenderState } from '../../render-state.js';
import { type Values, TextureSpec } from '../../renderable/schema.js';
import type { Texture } from '../../texture.js';
import { ShaderCode } from '../../shader-code.js';
import { ValueCell } from '@molstar/core/util';
import { unpackRGBToInt } from '@molstar/core/util/number-packing';
import { QuadSchema, QuadValues } from '../util.js';
import { isTimingMode } from '@molstar/core/util/debug';

export const HistopyramidSumSchema = {
  ...QuadSchema,
  tTexture: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, { sampleType: 'sint' }),
};
export const HistopyramidSumShaderCode = ShaderCode('sum');
type HistopyramidSumValues = Values<typeof HistopyramidSumSchema>;

const HistopyramidSumName = 'histopyramid-sum';

function getHistopyramidSumRenderable(ctx: WebGLContext, texture: Texture): ComputeRenderable<HistopyramidSumValues> {
  if (ctx.namedComputeRenderables[HistopyramidSumName]) {
    const v = ctx.namedComputeRenderables[HistopyramidSumName].values as HistopyramidSumValues;

    ValueCell.update(v.tTexture, texture);

    ctx.namedComputeRenderables[HistopyramidSumName].update();
  } else {
    ctx.namedComputeRenderables[HistopyramidSumName] = createHistopyramidSumRenderable(ctx, texture);
  }
  return ctx.namedComputeRenderables[HistopyramidSumName];
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

  if (!ctx.namedFramebuffers[HistopyramidSumName]) {
    ctx.namedFramebuffers[HistopyramidSumName] = resources.framebuffer();
  }
  const framebuffer = ctx.namedFramebuffers[HistopyramidSumName];

  if (!ctx.namedTextures[HistopyramidSumName]) {
    ctx.namedTextures[HistopyramidSumName] = capabilities.integerTextures
      ? resources.texture('image-int32', 'rgba', 'int', 'nearest')
      : resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    ctx.namedTextures[HistopyramidSumName].define(1, 1);
  }
  const sumTexture = ctx.namedTextures[HistopyramidSumName];
  sumTexture.attachFramebuffer(framebuffer, 0);

  ctx.beginRenderPass(framebuffer, { label: 'histogram pyramid sum' });
  ctx.setRenderState(RenderState.Default);
  ctx.disableScissor();
  ctx.setViewport(0, 0, 1, 1);
  renderable.render();
  gl.finish();

  ctx.readPixels(0, 0, 1, 1, capabilities.integerTextures ? sumInts : sumBytes);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  if (isTimingMode) ctx.timer.groupEnd('getHistopyramidSum');

  return capabilities.integerTextures ? sumInts[0] : unpackRGBToInt(sumBytes[0], sumBytes[1], sumBytes[2]);
}
