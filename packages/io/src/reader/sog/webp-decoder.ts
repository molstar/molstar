/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { SogImage } from './schema.js';

let gl: WebGL2RenderingContext | null | undefined;

function getContext() {
  if (gl === undefined) {
    // a context of its own, kept: images decode one after another
    gl = new OffscreenCanvas(1, 1).getContext('webgl2', { premultipliedAlpha: false, antialias: false });
  }
  return gl;
}

/**
 * Decodes the image `bytes` (lossless WebP) with the browser's decoder, the RGBA bytes as stored: decoded without
 * premultiplying or converting colors (`createImageBitmap`) and read back through a WebGL2 texture, not a 2D canvas,
 * which premultiplies (SOG's `sh0` keeps the opacity in alpha: its rgb, codebook indices, would change below 255).
 * Not in Node (no decoder)
 */
export async function decodeWebpInBrowser(bytes: Uint8Array, name: string): Promise<SogImage> {
  if (typeof createImageBitmap === 'undefined' || typeof OffscreenCanvas === 'undefined') {
    throw new Error(`decoding '${name}' needs a browser (WebP)`);
  }
  const bitmap = await createImageBitmap(new Blob([bytes as Uint8Array<ArrayBuffer>], { type: 'image/webp' }), {
    premultiplyAlpha: 'none',
    colorSpaceConversion: 'none',
  });
  const ctx = getContext();
  if (!ctx) {
    bitmap.close();
    throw new Error(`decoding '${name}' needs WebGL2`);
  }
  const { width, height } = bitmap;
  const texture = ctx.createTexture();
  const framebuffer = ctx.createFramebuffer();
  try {
    ctx.bindTexture(ctx.TEXTURE_2D, texture);
    ctx.pixelStorei(ctx.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
    ctx.pixelStorei(ctx.UNPACK_COLORSPACE_CONVERSION_WEBGL, ctx.NONE);
    ctx.pixelStorei(ctx.UNPACK_FLIP_Y_WEBGL, false);
    ctx.texImage2D(ctx.TEXTURE_2D, 0, ctx.RGBA8, ctx.RGBA, ctx.UNSIGNED_BYTE, bitmap);
    ctx.bindFramebuffer(ctx.FRAMEBUFFER, framebuffer);
    ctx.framebufferTexture2D(ctx.FRAMEBUFFER, ctx.COLOR_ATTACHMENT0, ctx.TEXTURE_2D, texture, 0);
    const data = new Uint8Array(width * height * 4);
    // the image's top row is the texture's first (no flip), read first
    ctx.readPixels(0, 0, width, height, ctx.RGBA, ctx.UNSIGNED_BYTE, data);
    return { width, height, data };
  } finally {
    ctx.bindFramebuffer(ctx.FRAMEBUFFER, null);
    ctx.deleteFramebuffer(framebuffer);
    ctx.deleteTexture(texture);
    bitmap.close();
  }
}
