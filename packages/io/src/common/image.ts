/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// An image's bytes (WebP, PNG, JPEG: any type the browser decodes) decoded by the browser as stored: not premultiplied,
// no color conversion. A format whose image channels are data, not colors (SOG's codebook indices, an opacity in
// alpha), needs exactly that: a 2D canvas premultiplies, changing the rgb of every pixel whose alpha is below 255; the
// readback goes through a WebGL texture instead (README.md §Decisions). Not in Node: it has no decoder.

/** A decoded image: RGBA bytes, row-major from the top-left, not premultiplied */
export interface DecodedImage {
  readonly width: number;
  readonly height: number;
  readonly data: Uint8Array;
}

/** An image's MIME type, as the browser's decoder takes it (`image/webp`, `image/png`) */
export type ImageType = `image/${string}`;

function assertBrowser(name: string) {
  if (typeof createImageBitmap === 'undefined' || typeof OffscreenCanvas === 'undefined') {
    throw new Error(`decoding the image '${name}' needs a browser`);
  }
}

/** The bitmap of `bytes` (of `type`) as stored: not premultiplied, no color conversion; for an upload as it is */
export function decodeImageBitmap(bytes: Uint8Array, type: ImageType, name = 'image'): Promise<ImageBitmap> {
  assertBrowser(name);
  return createImageBitmap(new Blob([bytes as Uint8Array<ArrayBuffer>], { type }), {
    premultiplyAlpha: 'none',
    colorSpaceConversion: 'none',
  });
}

let gl: WebGL2RenderingContext | WebGLRenderingContext | null | undefined;

/**
 * A context of its own, kept (images decode one after another): WebGL2, else WebGL1 (all the readback needs: unpacked
 * not premultiplied nor converted, an RGBA byte texture of any size without mipmaps a framebuffer's, `readPixels`)
 */
function getContext() {
  if (gl === undefined) {
    const canvas = new OffscreenCanvas(1, 1);
    const attributes = { premultipliedAlpha: false, antialias: false };
    gl = canvas.getContext('webgl2', attributes) ?? canvas.getContext('webgl', attributes);
  }
  return gl;
}

/**
 * The RGBA bytes of the image `bytes` (of `type`, the file `name`) as stored (`decodeImageBitmap`), read back through
 * a WebGL texture (a 2D canvas premultiplies)
 */
export async function decodeImage(
  bytes: Uint8Array,
  options: { readonly type: ImageType; readonly name: string },
): Promise<DecodedImage> {
  const { name } = options;
  const bitmap = await decodeImageBitmap(bytes, options.type, name);
  const ctx = getContext();
  if (!ctx) {
    bitmap.close();
    throw new Error(`decoding the image '${name}' needs WebGL`);
  }
  const { width, height } = bitmap;
  const texture = ctx.createTexture();
  const framebuffer = ctx.createFramebuffer();
  try {
    ctx.bindTexture(ctx.TEXTURE_2D, texture);
    ctx.pixelStorei(ctx.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
    ctx.pixelStorei(ctx.UNPACK_COLORSPACE_CONVERSION_WEBGL, ctx.NONE);
    ctx.pixelStorei(ctx.UNPACK_FLIP_Y_WEBGL, false);
    // unsized RGBA of bytes: WebGL1's (WebGL2's RGBA8 too)
    ctx.texImage2D(ctx.TEXTURE_2D, 0, ctx.RGBA, ctx.RGBA, ctx.UNSIGNED_BYTE, bitmap);
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
