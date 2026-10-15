/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { GraphicsContext } from '../context.js';
import type { WebGLContext } from '../webgl/context.js';
import { RenderState } from '../render-state.js';
import { createNullTexture, type Texture } from '../texture.js';
import { asGLTexture } from '../webgl/texture.js';
import { ValueCell } from '@molstar/core/util';
import { ValueSpec, AttributeSpec, UniformSpec, type Values, TextureSpec, DefineSpec } from '../renderable/schema.js';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { ShaderCode } from '../shader-code.js';
import { type ComputeRenderable, createComputeRenderable } from '../renderable.js';

export const QuadPositions = new Float32Array([
  1.0,
  1.0,
  -1.0,
  1.0,
  -1.0,
  -1.0, // First triangle
  -1.0,
  -1.0,
  1.0,
  -1.0,
  1.0,
  1.0, // Second triangle
]);

export const QuadSchema = {
  drawCount: ValueSpec('number'),
  instanceCount: ValueSpec('number'),
  aPosition: AttributeSpec('float32', 2, 0),
  uQuadScale: UniformSpec('v2'),
};

export const QuadValues: Values<typeof QuadSchema> = {
  drawCount: ValueCell.create(6),
  instanceCount: ValueCell.create(1),
  aPosition: ValueCell.create(QuadPositions),
  uQuadScale: ValueCell.create(Vec2.create(1, 1)),
};

//

export const CopySchema = {
  ...QuadSchema,
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),
  uTexSize: UniformSpec('v2'),
  /** pack the red channel (depth) into RGBA, for depth copies into non-float targets */
  dPackDepth: DefineSpec('boolean'),
};
export const CopyShaderCode = ShaderCode('copy');
export type CopyRenderable = ComputeRenderable<Values<typeof CopySchema>>;

export function createCopyRenderable(ctx: GraphicsContext, texture: Texture, packDepth = false): CopyRenderable {
  const values: Values<typeof CopySchema> = {
    ...QuadValues,
    tColor: ValueCell.create(texture),
    uTexSize: ValueCell.create(Vec2.create(texture.getWidth(), texture.getHeight())),
    dPackDepth: ValueCell.create(packDepth),
  };

  const schema = { ...CopySchema };
  const renderItem = ctx.createComputeRenderItem('triangles', CopyShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

const SharedCopyName = 'shared-copy';

export function getSharedCopyRenderable(ctx: GraphicsContext, texture: Texture) {
  if (!ctx.namedComputeRenderables[SharedCopyName]) {
    ctx.namedComputeRenderables[SharedCopyName] = createCopyRenderable(ctx, createNullTexture());
  }
  const copy = ctx.namedComputeRenderables[SharedCopyName] as CopyRenderable;
  ValueCell.update(copy.values.tColor, texture);
  ValueCell.update(
    copy.values.uTexSize,
    Vec2.set(copy.values.uTexSize.ref.value, texture.getWidth(), texture.getHeight()),
  );
  copy.update();
  return copy;
}

//

const ReadTextureName = 'read-texture';
const ReadSingleChannelTextureName = 'read-single-channel-texture';

export function readTexture<T extends Uint8Array | Float32Array | Int32Array = Uint8Array>(
  ctx: WebGLContext,
  texture: Texture,
  array?: T,
) {
  const { gl, resources } = ctx;
  if (!array && asGLTexture(texture).type !== gl.UNSIGNED_BYTE) throw new Error('unsupported texture type');

  if (!ctx.namedFramebuffers[ReadTextureName]) {
    ctx.namedFramebuffers[ReadTextureName] = resources.framebuffer();
  }
  const framebuffer = ctx.namedFramebuffers[ReadTextureName];

  const width = texture.getWidth();
  const height = texture.getHeight();
  if (!array) array = new Uint8Array(width * height * 4) as T;
  framebuffer.bind();
  texture.attachFramebuffer(framebuffer, 0);
  ctx.readPixels(0, 0, width, height, array);

  return { array, width, height };
}

/** Reads an `r` texture through a copy into `rgba` (single-channel formats aren't readable everywhere), values in `.r` */
export function readSingleChannelTexture(ctx: WebGLContext, texture: Texture) {
  const { gl, resources } = ctx;
  if (asGLTexture(texture).type !== gl.UNSIGNED_BYTE) throw new Error('unsupported texture type');

  const width = texture.getWidth();
  const height = texture.getHeight();

  const copy = getSharedCopyRenderable(ctx, texture);

  if (!ctx.namedFramebuffers[ReadSingleChannelTextureName]) {
    ctx.namedFramebuffers[ReadSingleChannelTextureName] = resources.framebuffer();
  }
  const framebuffer = ctx.namedFramebuffers[ReadSingleChannelTextureName];

  if (!ctx.namedTextures[ReadSingleChannelTextureName]) {
    ctx.namedTextures[ReadSingleChannelTextureName] = resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
  }
  const copyTex = ctx.namedTextures[ReadSingleChannelTextureName];
  copyTex.define(width, height);
  copyTex.attachFramebuffer(framebuffer, 0);

  ctx.beginRenderPass(framebuffer, { label: 'read single channel' });
  ctx.setRenderState(RenderState.Additive);
  ctx.setViewport(0, 0, width, height);
  ctx.setScissor(0, 0, width, height);
  ctx.clear({ color: [0, 0, 0, 0] });
  copy.render();

  const array = new Uint8Array(width * height * 4);
  ctx.readPixels(0, 0, width, height, array);

  return { array, width, height };
}
