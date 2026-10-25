/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { idFactory } from '@molstar/core/util/id-factory';
import type { ValueOf } from '@molstar/core/util/type-helpers';
import type { Framebuffer } from './framebuffer.js';
import type { GraphicsContext } from './context.js';
import type { DataArray, TextureImage, TextureVolume } from './renderable/util.js';

/** shared by all backends, texture ids are used as cache keys */
export const getNextTextureId = idFactory();

export type TextureKindValue = {
  'image-uint8': TextureImage<Uint8Array>;
  'image-float32': TextureImage<Float32Array>;
  'image-float16': TextureImage<Float32Array>;
  'image-int32': TextureImage<Int32Array>;
  'image-depth': TextureImage<Uint8Array>; // TODO should be Uint32Array
  'volume-uint8': TextureVolume<Uint8Array>;
  'volume-float32': TextureVolume<Float32Array>;
  'volume-float16': TextureVolume<Float32Array>;
  texture: Texture;
};
export type TextureValueType = ValueOf<TextureKindValue>;
export type TextureKind = keyof TextureKindValue;
/** for `depth`/`depth-stencil` formats a precision preference, backends fall back to the best supported type */
export type TextureType = 'ubyte' | 'ushort' | 'float' | 'fp16' | 'int' | 'float-stencil' | 'uint24-8';
/** `r` textures are sampled as `(value, 0, 0, 1)` or `(value, value, value, 1)`, read them from `.r` */
export type TextureFormat = 'r' | 'rg' | 'rgba' | 'depth' | 'depth-stencil';
/** Numbers are shortcuts for color attachment */
export type TextureAttachment =
  | 'depth'
  | 'stencil'
  | 'depth-stencil'
  | 'color0'
  | 'color1'
  | 'color2'
  | 'color3'
  | 'color4'
  | 'color5'
  | 'color6'
  | 'color7'
  | 0
  | 1
  | 2
  | 3
  | 4
  | 5
  | 6
  | 7;
export type TextureFilter = 'nearest' | 'linear';

export interface Texture {
  readonly id: number;

  getWidth: () => number;
  getHeight: () => number;
  getDepth: () => number;

  getByteCount: () => number;

  define: (width: number, height: number, depth?: number) => void;
  /**
   * The `sub` option requires an existing allocation on the GPU, that is, either
   * `define` or `load` without `sub` must have been called before.
   */
  load: (image: TextureImage<any> | TextureVolume<any> | HTMLImageElement, sub?: boolean) => void;
  /**
   * Writes the data array `data` (`DataSpec`) at element `offset` of the texture of `width` × `height` texels (allocated
   * at that size where it is another), an element per texel row by row: its element `ranges` (`[start, end)` pairs, see
   * `TextureChanges`), all where undefined. A backend's data texture (WebGL)
   */
  loadData?: (data: DataArray<any>, width: number, height: number, offset?: number, ranges?: readonly number[]) => void;
  /**
   * Writes the `rows` (`[start, end)` pairs) of the 2d `image` into the texture that holds an earlier state of it at its
   * size (`TextureChanges.writeRows`), the other rows kept; `image` stays the texture's source. Loads all of it where
   * the texture holds another size, or the image is flipped
   */
  loadRows?: (image: TextureImage<any>, rows: readonly number[]) => void;
  /** no-op where unsupported (non-power-of-two sizes on WebGL1) */
  mipmap: () => void;
  /** Use `layer` to attach a z-slice of a 3D texture */
  attachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment, layer?: number) => void;
  detachFramebuffer: (framebuffer: Framebuffer, attachment: TextureAttachment) => void;
  /**
   * Copies a region of the first color attachment of `framebuffer` to `dstX`, `dstY` of this 2D texture (rectangles
   * with a top-left origin), both must have the same format (a blit on WebGL2), not affected by scissor, write masks
   * or blending
   */
  copyFromFramebuffer: (
    framebuffer: Framebuffer,
    x: number,
    y: number,
    width: number,
    height: number,
    dstX: number,
    dstY: number,
  ) => void;

  /** How the texture was created and the data of its last `load`, see `transferTexture` */
  getSource?: () => TextureSource;

  /**
   * Frees the GPU storage, keeping the size: the next write allocates it again (`define`, lazy storage), the contents
   * are lost. For render targets of passes that stopped drawing into them.
   */
  releaseStorage: () => void;

  destroy: () => void;
}

/**
 * The creation arguments of a texture and the data of its last full `load` (undefined where the contents were rendered
 * or computed, or loaded in part)
 */
export interface TextureSource {
  readonly kind: TextureKind;
  readonly format: TextureFormat;
  readonly type: TextureType;
  readonly filter: TextureFilter;
  readonly data: TextureImage<any> | TextureVolume<any> | HTMLImageElement | undefined;
  readonly mipmap: boolean;
}

/**
 * A texture of `graphics` with the contents of `texture` (of another context, e.g. after a backend switch) where they
 * were loaded from data, undefined where they were rendered or computed
 */
export function transferTexture(texture: Texture, graphics: GraphicsContext): Texture | undefined {
  const source = texture.getSource?.();
  if (!source?.data) return;
  const t = graphics.resources.texture(source.kind, source.format, source.type, source.filter);
  t.load(source.data);
  if (source.mipmap) t.mipmap();
  return t;
}

export type TextureValues = { [k: string]: ValueCell<TextureValueType> };
export type Textures = [string, Texture][];

export function getBytesPerElement(format: TextureFormat, type: TextureType): number {
  return getFormatSize(format) * getTypeSize(type);
}

function getFormatSize(format: TextureFormat) {
  switch (format) {
    case 'r':
      return 1;
    case 'rg':
      return 2;
    case 'rgba':
      return 4;
    case 'depth':
      return 1;
    case 'depth-stencil':
      return 2;
  }
}

function getTypeSize(type: TextureType): number {
  switch (type) {
    case 'ubyte':
      return 1;
    case 'ushort':
      return 2;
    case 'float':
      return 4;
    case 'fp16':
      return 2;
    case 'int':
      return 4;
    case 'float-stencil':
      return 4;
    case 'uint24-8':
      return 2;
  }
}

//

const NullTextures = new WeakSet<Texture>();

export function isNullTexture(texture: Texture) {
  return NullTextures.has(texture);
}

/** Placeholder without GPU allocation, backends bind it as an empty texture */
export function createNullTexture(): Texture {
  const texture: Texture = {
    id: getNextTextureId(),

    getWidth: () => 0,
    getHeight: () => 0,
    getDepth: () => 0,
    getByteCount: () => 0,

    define: () => {},
    load: () => {},
    mipmap: () => {},
    attachFramebuffer: () => {
      throw new Error('cannot attach null-texture to a framebuffer');
    },
    detachFramebuffer: () => {
      throw new Error('cannot detach null-texture from a framebuffer');
    },
    copyFromFramebuffer: () => {
      throw new Error('cannot copy to null-texture');
    },

    releaseStorage: () => {},
    destroy: () => {},
  };
  NullTextures.add(texture);
  return texture;
}

/**
 * Decodes the image at `source` (a URL or data URL), see `GraphicsResources.loadImage`.
 * This will not work on node.js without a polyfill for `HTMLImageElement`.
 */
export async function decodeImage(source: string): Promise<HTMLImageElement> {
  const img = new Image();
  img.src = source;
  await img.decode();
  return img;
}

export type CubeSide = 'nx' | 'ny' | 'nz' | 'px' | 'py' | 'pz';

export type CubeFaces = {
  [k in CubeSide]: string | File | Promise<Blob>;
};
