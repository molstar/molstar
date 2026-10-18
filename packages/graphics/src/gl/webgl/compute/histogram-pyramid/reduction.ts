/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { WebGLContext } from '../../context.js';
import { createComputeRenderable, type ComputeRenderable } from '../../../renderable.js';
import { type Values, TextureSpec, UniformSpec } from '../../../renderable/schema.js';
import type { Texture, TextureFilter, TextureFormat, TextureKind, TextureType } from '../../../texture.js';
import { ShaderCode } from '../../../shader-code.js';
import { ValueCell } from '@molstar/core/util';
import { QuadSchema, QuadValues } from '../../../compute/util.js';
import type { Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { getHistopyramidSum } from './sum.js';
import type { Framebuffer } from '../../../framebuffer.js';
import { isPowerOfTwo } from '@molstar/core/math/misc';
import { RenderState } from '../../../render-state.js';
import { isTimingMode } from '@molstar/core/util/debug';
import { setTexelScissor, setTexelViewport } from '../../render-pass-target.js';

export const HistopyramidReductionSchema = {
  ...QuadSchema,
  tInputLevel: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tPreviousLevel: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, { sampleType: 'sint' }),
  uSize: UniformSpec('f'),
  uTexSize: UniformSpec('f'),
  uFirst: UniformSpec('b'),
};
export const HistopyramidReductionShaderCode = ShaderCode('reduction');
type HistopyramidReductionValues = Values<typeof HistopyramidReductionSchema>;

const HistogramPyramidName = 'histogram-pyramid';

function getHistopyramidReductionRenderable(
  ctx: WebGLContext,
  inputLevel: Texture,
  previousLevel: Texture,
): ComputeRenderable<HistopyramidReductionValues> {
  if (ctx.namedComputeRenderables[HistogramPyramidName]) {
    const v = ctx.namedComputeRenderables[HistogramPyramidName].values as HistopyramidReductionValues;

    ValueCell.update(v.tInputLevel, inputLevel);
    ValueCell.update(v.tPreviousLevel, previousLevel);

    ctx.namedComputeRenderables[HistogramPyramidName].update();
  } else {
    ctx.namedComputeRenderables[HistogramPyramidName] = createHistopyramidReductionRenderable(
      ctx,
      inputLevel,
      previousLevel,
    );
  }
  return ctx.namedComputeRenderables[HistogramPyramidName];
}

function createHistopyramidReductionRenderable(ctx: WebGLContext, inputLevel: Texture, previousLevel: Texture) {
  const values: HistopyramidReductionValues = {
    ...QuadValues,
    tInputLevel: ValueCell.create(inputLevel),
    tPreviousLevel: ValueCell.create(previousLevel),
    uSize: ValueCell.create(0),
    uTexSize: ValueCell.create(0),
    uFirst: ValueCell.create(true),
  };

  const schema = { ...HistopyramidReductionSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', HistopyramidReductionShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

type TextureFramebuffer = { texture: Texture; framebuffer: Framebuffer };
function getLevelTextureFramebuffer(ctx: WebGLContext, level: number) {
  const size = Math.pow(2, level);
  const name = `level${level}`;
  const texture = ctx.capabilities.integerTextures
    ? getTexture(name, ctx, 'image-int32', 'r', 'int', 'nearest')
    : getTexture(name, ctx, 'image-uint8', 'rgba', 'ubyte', 'nearest');
  texture.define(size, size);
  let framebuffer = tryGetFramebuffer(name, ctx);
  if (!framebuffer) {
    framebuffer = getFramebuffer(name, ctx);
  }
  texture.attachFramebuffer(framebuffer, 0);
  return { texture, framebuffer };
}

function getFramebuffer(name: string, webgl: WebGLContext): Framebuffer {
  const _name = `${HistogramPyramidName}-${name}`;
  if (!webgl.namedFramebuffers[_name]) {
    webgl.namedFramebuffers[_name] = webgl.resources.framebuffer();
  }
  return webgl.namedFramebuffers[_name];
}

function getTexture(
  name: string,
  webgl: WebGLContext,
  kind: TextureKind,
  format: TextureFormat,
  type: TextureType,
  filter: TextureFilter,
): Texture {
  const _name = `${HistogramPyramidName}-${name}`;
  if (!webgl.namedTextures[_name]) {
    webgl.namedTextures[_name] = webgl.resources.texture(kind, format, type, filter);
  }
  return webgl.namedTextures[_name];
}

function tryGetFramebuffer(name: string, webgl: WebGLContext): Framebuffer | undefined {
  const _name = `${HistogramPyramidName}-${name}`;
  return webgl.namedFramebuffers[_name];
}

export interface HistogramPyramid {
  pyramidTex: Texture;
  count: number;
  height: number;
  levels: number;
  scale: Vec2;
}

export function createHistogramPyramid(
  ctx: WebGLContext,
  inputTexture: Texture,
  scale: Vec2,
  gridTexDim: Vec3,
): HistogramPyramid {
  if (isTimingMode) ctx.timer.group('createHistogramPyramid');
  const w = inputTexture.getWidth();
  const h = inputTexture.getHeight();

  // printTexture(ctx, inputTexture, 2)
  if (w !== h || !isPowerOfTwo(w)) {
    throw new Error('inputTexture must be of square power-of-two size');
  }

  // This part set the levels
  const levels = Math.ceil(Math.log(w) / Math.log(2));
  const maxSize = Math.pow(2, levels);
  const maxSizeX = Math.pow(2, levels);
  const maxSizeY = Math.pow(2, levels - 1);
  // console.log('levels', levels, 'maxSize', maxSize, [maxSizeX, maxSizeY], 'input', w);

  const pyramidTex = ctx.capabilities.integerTextures
    ? getTexture('pyramid', ctx, 'image-int32', 'r', 'int', 'nearest')
    : getTexture('pyramid', ctx, 'image-uint8', 'rgba', 'ubyte', 'nearest');
  pyramidTex.define(maxSizeX, maxSizeY);

  const framebuffer = getFramebuffer('pyramid', ctx);
  pyramidTex.attachFramebuffer(framebuffer, 0);

  ctx.beginRenderPass(framebuffer, { label: 'histogram pyramid' });
  ctx.setRenderState(RenderState.Default);
  setTexelViewport(ctx, 0, 0, maxSizeX, maxSizeY);
  setTexelScissor(ctx, 0, 0, maxSizeX, maxSizeY);
  ctx.clear({ color: [0, 0, 0, 0] });

  const levelTexturesFramebuffers: TextureFramebuffer[] = [];
  for (let i = 0; i < levels; ++i) levelTexturesFramebuffers.push(getLevelTextureFramebuffer(ctx, i));

  const renderable = getHistopyramidReductionRenderable(ctx, inputTexture, levelTexturesFramebuffers[0].texture);

  let offset = 0;
  for (let i = 0; i < levels; i++) {
    const currLevel = levels - 1 - i;
    const tf = levelTexturesFramebuffers[currLevel];
    ctx.beginRenderPass(tf.framebuffer, { label: 'histogram pyramid' });

    const size = Math.pow(2, currLevel);
    // console.log('size', size, 'draw-level', currLevel, 'read-level', levels - i);

    ValueCell.update(renderable.values.uSize, Math.pow(2, i + 1) / maxSize);
    ValueCell.update(renderable.values.uTexSize, size);
    ValueCell.updateIfChanged(renderable.values.uFirst, i === 0);
    if (i > 0) {
      ValueCell.update(renderable.values.tPreviousLevel, levelTexturesFramebuffers[levels - i].texture);
      renderable.update();
    }
    setTexelViewport(ctx, 0, 0, size, size);
    setTexelScissor(ctx, 0, 0, size, size);
    ctx.clear({ color: [0, 0, 0, 0] });
    setTexelScissor(ctx, 0, 0, gridTexDim[0], gridTexDim[1]);
    renderable.render();

    // at texel row 0 (rectangles have a top-left origin)
    pyramidTex.copyFromFramebuffer(
      tf.framebuffer,
      0,
      tf.texture.getHeight() - size,
      size,
      size,
      offset,
      pyramidTex.getHeight() - size,
    );

    offset += size;
  }

  if (isTimingMode) ctx.timer.groupEnd('createHistogramPyramid');

  // printTextureImage(readTexture(ctx, pyramidTex), { scale: 0.75 });

  //

  // return at least a count of one to avoid issues downstram
  const count = Math.max(1, getHistopyramidSum(ctx, levelTexturesFramebuffers[0].texture));
  const height = Math.ceil(count / Math.pow(2, levels));
  // console.log({ height, count, scale });

  return { pyramidTex, count, height, levels, scale };
}
