/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SegmentedArray } from '../segmented-array.js';
import type { WebGLContext } from './context.js';
import type { DataArray, TextureImage, TextureVolume } from '../renderable/util.js';
import type { RenderableSchema } from '../renderable/schema.js';
import type { SharedTexture, SharedTextureImage } from '../shared-texture.js';
import type { Framebuffer } from '../framebuffer.js';
import {
  bindForAttachment,
  getFramebufferHeight,
  setAttachmentStorage,
  setFramebufferHeight,
  setIntegerColorAttachment,
} from './framebuffer.js';
import { endRenderPass } from './render-pass.js';
import type { WebGLState } from './state.js';
import { isWebGL2, type GLRenderingContext } from './compat.js';
import { isPromiseLike } from '@molstar/core/util/type-helpers';
import type { WebGLExtensions } from './extensions.js';
import { objectForEach } from '@molstar/core/util/object';
import { isPowerOfTwo } from '@molstar/core/math/misc';
import {
  type CubeFaces,
  type CubeSide,
  getBytesPerElement,
  getNextTextureId,
  type Texture,
  type TextureAttachment,
  type TextureFilter,
  type TextureFormat,
  type TextureKind,
  type Textures,
  type TextureType,
  type TextureValues,
} from '../texture.js';

export function getTarget(gl: GLRenderingContext, kind: TextureKind): number {
  switch (kind) {
    case 'image-uint8':
      return gl.TEXTURE_2D;
    case 'image-float32':
      return gl.TEXTURE_2D;
    case 'image-float16':
      return gl.TEXTURE_2D;
    case 'image-depth':
      return gl.TEXTURE_2D;
  }
  if (isWebGL2(gl)) {
    switch (kind) {
      case 'image-int32':
        return gl.TEXTURE_2D;
      case 'volume-uint8':
        return gl.TEXTURE_3D;
      case 'volume-float32':
        return gl.TEXTURE_3D;
      case 'volume-float16':
        return gl.TEXTURE_3D;
    }
  }
  throw new Error(`unknown texture kind '${kind}'`);
}

export function getFormat(gl: GLRenderingContext, format: TextureFormat, type: TextureType): number {
  switch (format) {
    case 'r':
      // WebGL1 has no red formats, luminance samples as `(L, L, L, 1)`
      if (!isWebGL2(gl)) return gl.LUMINANCE;
      return type === 'int' ? gl.RED_INTEGER : gl.RED;
    case 'rg':
      if (isWebGL2(gl) && (type === 'float' || type === 'fp16')) return gl.RG;
      else if (isWebGL2(gl) && type === 'int') return gl.RG_INTEGER;
      else throw new Error('texture format "rg" requires webgl2 and type "float" or int"');
    case 'rgba':
      if (isWebGL2(gl) && type === 'int') return gl.RGBA_INTEGER;
      return gl.RGBA;
    case 'depth':
      return gl.DEPTH_COMPONENT;
    case 'depth-stencil':
      return (gl as WebGL2RenderingContext).DEPTH_STENCIL;
  }
}

export function getInternalFormat(gl: GLRenderingContext, format: TextureFormat, type: TextureType): number {
  if (isWebGL2(gl)) {
    switch (format) {
      case 'r':
        switch (type) {
          case 'ubyte':
            return gl.R8;
          case 'float':
            return gl.R32F;
          case 'fp16':
            return gl.R16F;
          case 'int':
            return gl.R32I;
        }
      case 'rg':
        switch (type) {
          case 'ubyte':
            return gl.RG;
          case 'float':
            return gl.RG32F;
          case 'fp16':
            return gl.RG16F;
          case 'int':
            return gl.RG32I;
        }
      case 'rgba':
        switch (type) {
          case 'ubyte':
            return gl.RGBA;
          case 'float':
            return gl.RGBA32F;
          case 'fp16':
            return gl.RGBA16F;
          case 'int':
            return gl.RGBA32I;
        }
      case 'depth':
        switch (type) {
          case 'ushort':
            return gl.DEPTH_COMPONENT16;
          case 'float':
            return gl.DEPTH_COMPONENT32F;
        }
        break;
      case 'depth-stencil':
        switch (type) {
          case 'float-stencil':
            return gl.DEPTH32F_STENCIL8;
          case 'uint24-8':
            return gl.DEPTH24_STENCIL8;
        }
        break;
    }
  }
  return getFormat(gl, format, type);
}

function getByteCount(format: TextureFormat, type: TextureType, width: number, height: number, depth: number): number {
  return getBytesPerElement(format, type) * width * height * (depth || 1);
}

export function getType(gl: GLRenderingContext, extensions: WebGLExtensions, type: TextureType): number {
  switch (type) {
    case 'ubyte':
      return gl.UNSIGNED_BYTE;
    case 'ushort':
      return gl.UNSIGNED_SHORT;
    case 'float':
      return gl.FLOAT;
    case 'fp16':
      if (extensions.textureHalfFloat) return extensions.textureHalfFloat.HALF_FLOAT;
      else throw new Error('extension "texture_half_float" unavailable');
    case 'int':
      if (isWebGL2(gl)) return gl.INT;
      else throw new Error('texture type "int" requires webgl2');
    case 'float-stencil':
      if (isWebGL2(gl)) return gl.FLOAT_32_UNSIGNED_INT_24_8_REV;
      else throw new Error('texture type "float-stencil" requires webgl2');
    case 'uint24-8':
      if (extensions.depthTexture) return extensions.depthTexture.UNSIGNED_INT_24_8;
      else throw new Error('extension "depth_texture" unavailable');
  }
}

export function getFilter(gl: GLRenderingContext, type: TextureFilter): number {
  switch (type) {
    case 'nearest':
      return gl.NEAREST;
    case 'linear':
      return gl.LINEAR;
  }
}

export function getAttachment(
  gl: GLRenderingContext,
  extensions: WebGLExtensions,
  attachment: TextureAttachment,
): number {
  switch (attachment) {
    case 'depth':
      return gl.DEPTH_ATTACHMENT;
    case 'stencil':
      return gl.STENCIL_ATTACHMENT;
    case 'depth-stencil':
      return (gl as WebGL2RenderingContext).DEPTH_STENCIL_ATTACHMENT;
    case 'color0':
    case 0:
      return gl.COLOR_ATTACHMENT0;
  }
  if (extensions.drawBuffers) {
    switch (attachment) {
      case 'color1':
      case 1:
        return extensions.drawBuffers.COLOR_ATTACHMENT1;
      case 'color2':
      case 2:
        return extensions.drawBuffers.COLOR_ATTACHMENT2;
      case 'color3':
      case 3:
        return extensions.drawBuffers.COLOR_ATTACHMENT3;
      case 'color4':
      case 4:
        return extensions.drawBuffers.COLOR_ATTACHMENT4;
      case 'color5':
      case 5:
        return extensions.drawBuffers.COLOR_ATTACHMENT5;
      case 'color6':
      case 6:
        return extensions.drawBuffers.COLOR_ATTACHMENT6;
      case 'color7':
      case 7:
        return extensions.drawBuffers.COLOR_ATTACHMENT7;
    }
  }
  throw new Error('unknown texture attachment');
}

function isImage(x: TextureImage<any> | TextureVolume<any> | HTMLImageElement): x is HTMLImageElement {
  return typeof HTMLImageElement !== 'undefined' && x instanceof HTMLImageElement;
}

/** Bytes of rows composed per upload of a segmented data array */
const SegmentedChunkBytes = 4 * 1024 * 1024;

function isTexture2d(
  x: TextureImage<any> | TextureVolume<any>,
  target: number,
  gl: GLRenderingContext,
): x is TextureImage<any> {
  return target === gl.TEXTURE_2D;
}

function isTexture3d(
  x: TextureImage<any> | TextureVolume<any>,
  target: number,
  gl: WebGL2RenderingContext,
): x is TextureVolume<any> {
  return target === gl.TEXTURE_3D;
}

/** WebGL texture, all textures created by the WebGL backend implement it */
export interface GLTexture extends Texture {
  readonly target: number;
  readonly format: number;
  readonly internalFormat: number;
  readonly type: number;
  readonly filter: number;

  bind: (id: TextureId) => void;
  unbind: (id: TextureId) => void;
}

export type TextureId = 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15;

/** Textures handed out by the WebGL resources, see `isGLTexture` */
const GLTextures = new WeakSet<Texture>();

/** Marks `texture` as created by the WebGL backend, called on the object the resources return */
export function registerGLTexture<T extends GLTexture>(texture: T): T {
  GLTextures.add(texture);
  return texture;
}

export function isGLTexture(texture: Texture): texture is GLTexture {
  return GLTextures.has(texture);
}

/** Narrows to `GLTexture`, throws for textures of another backend */
export function asGLTexture(texture: Texture): GLTexture {
  if (!isGLTexture(texture)) throw new Error('expected a WebGL texture, got a texture of another backend');
  return texture;
}

function getTexture(gl: GLRenderingContext) {
  const texture = gl.createTexture();
  if (texture === null) {
    throw new Error('Could not create WebGL texture');
  }
  return texture;
}

export function createTexture(
  gl: GLRenderingContext,
  extensions: WebGLExtensions,
  kind: TextureKind,
  _format: TextureFormat,
  _type: TextureType,
  _filter: TextureFilter,
  state: WebGLState,
): GLTexture {
  const id = getNextTextureId();
  const texture = getTexture(gl);

  // depth types are precision preferences, WebGL1 falls back to what `WEBGL_depth_texture` supports
  if (!isWebGL2(gl) && (_format === 'depth' || _format === 'depth-stencil')) {
    if (_type === 'float') _type = 'ushort';
    else if (_type === 'float-stencil') _type = 'uint24-8';
  }

  // check texture kind and type compatability
  if (
    (kind.endsWith('float32') && _type !== 'float') ||
    (kind.endsWith('float16') && _type !== 'fp16') ||
    (kind.endsWith('uint8') && _type !== 'ubyte') ||
    (kind.endsWith('int32') && _type !== 'int') ||
    (kind.endsWith('depth') &&
      _type !== 'ushort' &&
      _type !== 'float' &&
      _type !== 'float-stencil' &&
      _type !== 'uint24-8')
  ) {
    throw new Error(`texture kind '${kind}' and type '${_type}' are incompatible`);
  }

  if (!extensions.depthTexture && (_format === 'depth' || _format === 'depth-stencil')) {
    throw new Error(`extension 'WEBGL_depth_texture' needed for 'depth' texture format`);
  }

  const target = getTarget(gl, kind);
  const filter = getFilter(gl, _filter);
  const format = getFormat(gl, _format, _type);
  const internalFormat = getInternalFormat(gl, _format, _type);
  const type = getType(gl, extensions, _type);

  /** binding a texture replaces the one on the active unit, render items have to rebind theirs */
  function invalidateBindings() {
    state.currentRenderItemId = -1;
    state.boundTextures.length = 0;
  }

  function init() {
    invalidateBindings();
    gl.bindTexture(target, texture);
    gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, filter);
    gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, filter);
    // clamp-to-edge needed for non-power-of-two textures in webgl
    gl.texParameteri(target, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(target, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.bindTexture(target, null);
  }
  init();

  let width = 0,
    height = 0,
    depth = 0;
  let hasMipmap = false;
  /** the data of the last load covering the contents (not after a partial load or a later `define` of another size) */
  let sourceData: undefined | TextureImage<any> | TextureVolume<any> | HTMLImageElement;
  let destroyed = false;
  /** draw framebuffer with this texture attached, for `copyFromFramebuffer` on WebGL2 */
  let blitFramebuffer: WebGLFramebuffer | null = null;

  /** the size of the GL image (the storage), the defined size once allocated, see `ensureStorage` */
  let storedWidth = 0,
    storedHeight = 0,
    storedDepth = 0;
  const is3d = isWebGL2(gl) && target === gl.TEXTURE_3D;

  /** Specifies the GL image of the given size, contents undefined */
  function specify(w: number, h: number, d: number) {
    invalidateBindings();
    gl.bindTexture(target, texture);
    if (target === gl.TEXTURE_2D) {
      gl.texImage2D(target, 0, internalFormat, w, h, 0, format, type, null);
    } else if (is3d) {
      gl.texImage3D(target, 0, internalFormat, w, h, d, 0, format, type, null);
    } else {
      throw new Error('unknown texture target');
    }
    // not left bound at the active unit: a sampler no pass provides would read it, also while it is drawn into
    gl.bindTexture(target, null);
    (storedWidth = w), (storedHeight = h), (storedDepth = d);
  }

  /**
   * Records the size; the storage is allocated when the texture is first used as an attachment of a bound framebuffer,
   * copied into or partially loaded (`ensureStorage`). Until then the GL image is 1×1: sampled, it reads zero as an
   * undefined texture.
   */
  function define(_width: number, _height: number, _depth?: number) {
    if (_width === 0 || _height === 0 || (is3d && _depth === 0)) {
      throw new Error('empty textures are not allowed');
    }

    if (width === _width && height === _height && depth === (_depth || 0)) return;

    sourceData = undefined;
    (width = _width), (height = _height), (depth = _depth || 0);
    // storage of another size released
    const d = is3d ? 1 : 0;
    if (storedWidth !== 1 || storedHeight !== 1 || storedDepth !== d) specify(1, 1, d);
  }
  define(1, 1, is3d ? 1 : 0);

  /** Allocates the storage of the defined size, see `define` */
  function ensureStorage() {
    if (destroyed) return;
    if (storedWidth !== width || storedHeight !== height || storedDepth !== depth) specify(width, height, depth);
  }

  /** components per texel of a data array: a data element per texel */
  const dataItemSize = _format === 'r' ? 1 : 4;

  /**
   * Writes `data` (`DataSpec`) at element `offset` of the texture of `w` × `h` texels (allocated at that size where it is
   * another), an element per texel row by row: the element `ranges` of `data` (`[start, end)` pairs), all where
   * undefined. Whole rows from views of the array, the partial rows at the ends apart, no padded copy (composed from the
   * parts where segmented, a few rows at a time)
   */
  function loadData(data: DataArray<any>, w: number, h: number, offset = 0, ranges?: readonly number[]) {
    const array = data.array as (ArrayBufferView & ArrayLike<number>) | SegmentedArray;
    const count = data.count;
    invalidateBindings();
    if (storedWidth !== w || storedHeight !== h) {
      (width = w), (height = h), (depth = 0);
      specify(w, h, 0);
    }
    gl.bindTexture(target, texture);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    const segmented = SegmentedArray.is(array);
    const read = (start: number, end: number) =>
      segmented
        ? SegmentedArray.read(array, start * dataItemSize, end * dataItemSize)
        : (array as unknown as Uint8Array).subarray(start * dataItemSize, end * dataItemSize);
    const chunkRows = segmented
      ? Math.max(1, Math.floor(SegmentedChunkBytes / (w * dataItemSize * array.BYTES_PER_ELEMENT)))
      : h;
    const upload = (start: number, end: number) => {
      end = Math.min(end, count);
      let s = start;
      if (end <= s) return;
      // the elements before the first whole row
      const x = (offset + s) % w;
      if (x !== 0) {
        const n = Math.min(end - s, w - x);
        gl.texSubImage2D(target, 0, x, Math.floor((offset + s) / w), n, 1, format, type, read(s, s + n));
        s += n;
      }
      while (end - s >= w) {
        const rows = Math.min(chunkRows, Math.floor((end - s) / w));
        gl.texSubImage2D(target, 0, 0, (offset + s) / w, w, rows, format, type, read(s, s + rows * w));
        s += rows * w;
      }
      if (end > s) gl.texSubImage2D(target, 0, 0, (offset + s) / w, end - s, 1, format, type, read(s, end));
    };
    if (ranges) for (let k = 0; k < ranges.length; k += 2) upload(ranges[k], ranges[k + 1]);
    else upload(0, count);
    gl.bindTexture(target, null);
    sourceData = undefined;
  }

  function load(data: TextureImage<any> | TextureVolume<any> | HTMLImageElement, sub = false) {
    if (
      data.width === 0 ||
      data.height === 0 ||
      (!isImage(data) && isWebGL2(gl) && isTexture3d(data, target, gl) && data.depth === 0)
    ) {
      throw new Error('empty textures are not allowed');
    }

    // partial loads write into the storage of the defined size
    if (sub) ensureStorage();
    invalidateBindings();
    gl.bindTexture(target, texture);
    // unpack alignment of 1 since we use textures only for data
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
    gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, 0);
    if (isImage(data)) {
      (width = data.width), (height = data.height);
      (storedWidth = width), (storedHeight = height);
      gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
      gl.bindTexture(gl.TEXTURE_2D, texture);
      gl.texImage2D(gl.TEXTURE_2D, 0, internalFormat, format, type, data);
    } else if (isTexture2d(data, target, gl)) {
      const _filter = data.filter ? getFilter(gl, data.filter) : filter;
      gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, _filter);
      gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, _filter);
      gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, !!data.flipY);
      if (sub) {
        gl.texSubImage2D(target, 0, 0, 0, data.width, data.height, format, type, data.array);
      } else {
        (width = data.width), (height = data.height);
        (storedWidth = width), (storedHeight = height), (storedDepth = depth);
        gl.texImage2D(target, 0, internalFormat, width, height, 0, format, type, data.array);
      }
    } else if (isWebGL2(gl) && isTexture3d(data, target, gl)) {
      gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
      if (sub) {
        gl.texSubImage3D(target, 0, 0, 0, 0, data.width, data.height, data.depth, format, type, data.array);
      } else {
        (width = data.width), (height = data.height), (depth = data.depth);
        (storedWidth = width), (storedHeight = height), (storedDepth = depth);
        gl.texImage3D(target, 0, internalFormat, width, height, depth, 0, format, type, data.array);
      }
    } else {
      throw new Error('unknown texture target');
    }
    gl.bindTexture(target, null);
    sourceData = sub ? undefined : data;
  }

  /** `Texture.loadRows`: the rows from views of the array, as `loadData` writes */
  function loadRows(image: TextureImage<any>, rows: readonly number[]) {
    if (target !== gl.TEXTURE_2D || image.flipY || storedWidth !== image.width || storedHeight !== image.height) {
      load(image);
      return;
    }
    const array = image.array as Uint8Array;
    const stride = array.length / image.height;
    invalidateBindings();
    gl.bindTexture(target, texture);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    for (let k = 0; k < rows.length; k += 2) {
      const start = rows[k];
      const end = Math.min(rows[k + 1], image.height);
      if (end <= start) continue;
      const data = array.subarray(start * stride, end * stride);
      gl.texSubImage2D(target, 0, 0, start, image.width, end - start, format, type, data);
    }
    gl.bindTexture(target, null);
    sourceData = image;
  }

  function mipmap() {
    if (target !== gl.TEXTURE_2D) {
      throw new Error('mipmap only supported for 2d textures');
    }

    if (isWebGL2(gl) || (isPowerOfTwo(width) && isPowerOfTwo(height))) {
      ensureStorage();
      invalidateBindings();
      gl.bindTexture(target, texture);
      gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
      gl.generateMipmap(target);
      gl.bindTexture(target, null);
      hasMipmap = true;
    }
  }

  function resolveAttachment(attachment: TextureAttachment): TextureAttachment {
    return _format === 'depth-stencil' && attachment === 'depth' ? 'depth-stencil' : attachment;
  }

  function attachFramebuffer(framebuffer: Framebuffer, attachment: TextureAttachment, layer?: number) {
    const att = resolveAttachment(attachment);
    endRenderPass(state);
    bindForAttachment(gl, framebuffer);
    setAttachmentStorage(gl, framebuffer, getAttachment(gl, extensions, att), { kind: 'texture', id }, ensureStorage);
    if (att === 'depth-stencil') {
      gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, null);
    }
    if (target === gl.TEXTURE_2D) {
      gl.framebufferTexture2D(gl.FRAMEBUFFER, getAttachment(gl, extensions, att), gl.TEXTURE_2D, texture, 0);
    } else if (isWebGL2(gl) && target === gl.TEXTURE_3D) {
      if (layer === undefined) throw new Error('need `layer` to attach 3D texture');
      gl.framebufferTextureLayer(gl.FRAMEBUFFER, getAttachment(gl, extensions, att), texture, 0, layer);
    } else {
      throw new Error('unknown/unsupported texture target');
    }
    if (att === 0 || att === 'color0') setIntegerColorAttachment(gl, _type === 'int');
    setFramebufferHeight(framebuffer, () => height);
  }

  return {
    id,
    target,
    format,
    internalFormat,
    type,
    filter,

    getWidth: () => width,
    getHeight: () => height,
    getDepth: () => depth,

    getByteCount: () => getByteCount(_format, _type, storedWidth, storedHeight, storedDepth),

    define,
    load,
    loadData,
    loadRows,
    mipmap,
    getSource: () => ({ kind, format: _format, type: _type, filter: _filter, data: sourceData, mipmap: hasMipmap }),
    bind: (id: TextureId) => {
      gl.activeTexture(gl.TEXTURE0 + id);
      gl.bindTexture(target, texture);
    },
    unbind: (id: TextureId) => {
      invalidateBindings();
      gl.activeTexture(gl.TEXTURE0 + id);
      gl.bindTexture(target, null);
    },
    attachFramebuffer,
    detachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => {
      const att = resolveAttachment(attachment);
      endRenderPass(state);
      bindForAttachment(gl, framebuffer);
      setAttachmentStorage(gl, framebuffer, getAttachment(gl, extensions, att), undefined);
      if (target === gl.TEXTURE_2D) {
        gl.framebufferTexture2D(gl.FRAMEBUFFER, getAttachment(gl, extensions, att), gl.TEXTURE_2D, null, 0);
      } else if (isWebGL2(gl) && target === gl.TEXTURE_3D) {
        gl.framebufferTextureLayer(gl.FRAMEBUFFER, getAttachment(gl, extensions, att), null, 0, 0);
      } else {
        throw new Error('unknown texture target');
      }
      if (att === 0 || att === 'color0') setIntegerColorAttachment(gl, false);
    },
    copyFromFramebuffer: (
      framebuffer: Framebuffer,
      x: number,
      srcY: number,
      w: number,
      h: number,
      dstX: number,
      dstTopY: number,
    ) => {
      if (target !== gl.TEXTURE_2D) throw new Error('copy only supported for 2d textures');
      // rectangles have a top-left origin, GL's a bottom-left one
      const y = getFramebufferHeight(framebuffer) - srcY - h;
      const dstY = height - dstTopY - h;
      endRenderPass(state);
      ensureStorage();
      framebuffer.bind();
      if (isWebGL2(gl)) {
        // blit into a framebuffer with this texture attached, leaves the texture units alone
        if (!blitFramebuffer) {
          blitFramebuffer = gl.createFramebuffer();
          gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, blitFramebuffer);
          gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, texture, 0);
        } else {
          gl.bindFramebuffer(gl.DRAW_FRAMEBUFFER, blitFramebuffer);
        }
        // blits are clipped by the scissor test, but not affected by write masks or blending
        const scissorTest = gl.isEnabled(gl.SCISSOR_TEST);
        if (scissorTest) gl.disable(gl.SCISSOR_TEST);
        gl.blitFramebuffer(x, y, x + w, y + h, dstX, dstY, dstX + w, dstY + h, gl.COLOR_BUFFER_BIT, gl.NEAREST);
        if (scissorTest) gl.enable(gl.SCISSOR_TEST);
        framebuffer.bind();
      } else {
        invalidateBindings();
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(target, texture);
        gl.copyTexSubImage2D(target, 0, dstX, dstY, x, y, w, h);
        gl.bindTexture(target, null);
      }
    },
    releaseStorage: () => {
      if (destroyed) return;
      sourceData = undefined;
      hasMipmap = false;
      // a 1×1 image until written again at the defined size (`ensureStorage`)
      const d = is3d ? 1 : 0;
      if (storedWidth !== 1 || storedHeight !== 1 || storedDepth !== d) specify(1, 1, d);
    },
    destroy: () => {
      if (destroyed) return;
      gl.deleteTexture(texture);
      if (blitFramebuffer) gl.deleteFramebuffer(blitFramebuffer);
      destroyed = true;
    },
  };
}

/** `shared` receives the textures of specs with `shared`, by name (`.v6/plans/gl/shared-textures.md` §4) */
export function createTextures(
  ctx: WebGLContext,
  schema: RenderableSchema,
  values: TextureValues,
  shared: Map<string, SharedTexture>,
) {
  const { resources } = ctx;
  const textures: Textures = [];
  Object.keys(schema).forEach((k) => {
    const spec = schema[k];
    if (spec.type === 'texture') {
      const value = values[k];
      if (value) {
        if (spec.kind === 'texture') {
          textures[textures.length] = [k, value.ref.value as Texture];
        } else if (spec.shared) {
          const held = resources.sharedTextures.acquire(spec, value.ref.value as SharedTextureImage);
          shared.set(k, held);
          textures[textures.length] = [k, held.texture];
        } else if (spec.data !== undefined) {
          // a data array (`DataSpec`): in a range of the data pages, which the render item binds
        } else {
          const texture = resources.texture(spec.kind, spec.format, spec.dataType, spec.filter);
          texture.load(value.ref.value as TextureImage<any> | TextureVolume<any>);
          textures[textures.length] = [k, texture];
        }
      }
    }
  });
  return textures;
}

//

export function getCubeTarget(gl: GLRenderingContext, side: CubeSide): number {
  switch (side) {
    case 'nx':
      return gl.TEXTURE_CUBE_MAP_NEGATIVE_X;
    case 'ny':
      return gl.TEXTURE_CUBE_MAP_NEGATIVE_Y;
    case 'nz':
      return gl.TEXTURE_CUBE_MAP_NEGATIVE_Z;
    case 'px':
      return gl.TEXTURE_CUBE_MAP_POSITIVE_X;
    case 'py':
      return gl.TEXTURE_CUBE_MAP_POSITIVE_Y;
    case 'pz':
      return gl.TEXTURE_CUBE_MAP_POSITIVE_Z;
  }
}

export function createCubeTexture(
  gl: GLRenderingContext,
  faces: CubeFaces,
  mipmaps: boolean,
  onload?: (errored?: boolean) => void,
): GLTexture {
  const target = gl.TEXTURE_CUBE_MAP;
  const filter = gl.LINEAR;
  const internalFormat = gl.RGBA;
  const format = gl.RGBA;
  const type = gl.UNSIGNED_BYTE;

  let size = 0;

  const texture = gl.createTexture();
  gl.bindTexture(target, texture);

  function load(cubeTarget: number, level: number, image: HTMLImageElement, isReset: boolean) {
    if (size === 0) size = image.width;

    gl.bindTexture(target, texture);
    gl.texImage2D(cubeTarget, level, internalFormat, size, size, 0, format, type, null);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);
    gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);
    gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, 0);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    gl.bindTexture(target, texture);
    gl.texImage2D(cubeTarget, level, internalFormat, format, type, image);

    loadedCount += 1;
    if (loadedCount === 6) {
      if (!destroyed) {
        if (mipmaps) {
          gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
          gl.generateMipmap(target);
        } else {
          gl.texParameteri(target, gl.TEXTURE_MIN_FILTER, filter);
        }
        gl.texParameteri(target, gl.TEXTURE_MAG_FILTER, filter);
      }
      if (!isReset) onload?.(destroyed);
    }
  }

  const facesData: { cubeTarget: number; level: number; image: HTMLImageElement }[] = [];

  let loadedCount = 0;
  objectForEach(faces, (source, side) => {
    if (!source) return;

    const level = 0;
    const cubeTarget = getCubeTarget(gl, side as CubeSide);

    const image = new Image();
    if (source instanceof File) {
      image.src = URL.createObjectURL(source);
    } else if (isPromiseLike(source)) {
      source.then((blob) => {
        image.src = URL.createObjectURL(blob);
      });
    } else {
      image.src = source;
    }

    facesData.push({ cubeTarget, level, image });

    image.addEventListener('load', () => {
      load(cubeTarget, level, image, false);
    });
    image.addEventListener('error', () => {
      onload?.(true);
    });
  });

  let destroyed = false;

  return {
    id: getNextTextureId(),
    target,
    format,
    internalFormat,
    type,
    filter,

    getWidth: () => size,
    getHeight: () => size,
    getDepth: () => 0,
    getByteCount: () => {
      return getByteCount('rgba', 'ubyte', size, size, 0) * 6 * (mipmaps ? 2 : 1);
    },

    define: () => {},
    load: () => {},
    mipmap: () => {},
    bind: (id: TextureId) => {
      gl.activeTexture(gl.TEXTURE0 + id);
      gl.bindTexture(target, texture);
    },
    unbind: (id: TextureId) => {
      gl.activeTexture(gl.TEXTURE0 + id);
      gl.bindTexture(target, null);
    },
    attachFramebuffer: () => {},
    detachFramebuffer: () => {},
    copyFromFramebuffer: () => {
      throw new Error('copy not supported for cube textures');
    },
    // loaded from images, not a target
    releaseStorage: () => {},

    destroy: () => {
      if (destroyed) return;
      gl.deleteTexture(texture);
      destroyed = true;
    },
  };
}
