/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { idFactory } from '@molstar/core/util/id-factory';
import type { TextureImage, TextureVolume } from '../renderable/util.js';
import type { Texture, TextureFilter, TextureFormat, TextureKind, TextureType } from '../texture.js';
import { TextureUsage } from './constants.js';
import { notSupported } from './errors.js';
import { asWebGPUFramebuffer } from './framebuffer.js';
import { getMipLevelCount, type MipmapGenerator } from './mipmap.js';

const getNextTextureId = idFactory();

/**
 * Destroys a GPU resource once no recorded commands use it: WebGPU invalidates destroyed resources for
 * all later submits, also for commands recorded before
 */
export type Release = (resource: { destroy(): void }) => void;

/** A resource recreated on the current device by `reset` when the context restores a lost device */
export interface WebGPURestorable {
  reset: () => void;
}

/** What WebGPU resources need from their context */
export interface WebGPUResourceEnv {
  /** the current device, read it on use: it changes when the context restores a lost device */
  readonly device: GPUDevice;
  /** resources reset after a device was restored, resources remove themselves when destroyed */
  readonly restorables: Set<WebGPURestorable>;
  readonly release: Release;
  /** Encodes commands outside of render passes (copies), ends the open pass */
  readonly encode: (commands: (encoder: GPUCommandEncoder) => void) => void;
  readonly generateMipmaps: MipmapGenerator;
}

export interface WebGPUTexture extends Texture {
  readonly format: GPUTextureFormat;
  readonly dimension: GPUTextureDimension;
  /** view for sampling (depth aspect of depth-stencil formats), undefined before `define` or `load` */
  getView: () => GPUTextureView | undefined;
  /** view for render pass attachments (all aspects), undefined before `define` or `load` */
  getAttachmentView: () => GPUTextureView | undefined;
  getGPUTexture: () => GPUTexture | undefined;
  readonly hasStencil: boolean;
  /** as on WebGL, a property of the texture; programs sample linear textures filtered, see `isLinearFiltered` */
  readonly filter: TextureFilter;
  /** incremented when the GPU texture is recreated, bind groups referencing the old view are stale */
  readonly version: number;
  /** as `Texture.load`, also from image bitmaps (decoded without color space conversion, see `loadImage`) */
  load: (image: TextureImage<any> | TextureVolume<any> | HTMLImageElement | ImageBitmap, sub?: boolean) => void;
}

const WebGPUTextures = new WeakSet<Texture>();

export function isWebGPUTexture(texture: Texture): texture is WebGPUTexture {
  return WebGPUTextures.has(texture);
}

/** Narrows to `WebGPUTexture`, throws for textures of another backend */
export function asWebGPUTexture(texture: Texture): WebGPUTexture {
  if (!isWebGPUTexture(texture)) throw new Error('expected a WebGPU texture, got a texture of another backend');
  return texture;
}

export function getGPUTextureFormat(format: TextureFormat, type: TextureType): GPUTextureFormat {
  switch (format) {
    case 'depth':
      return type === 'float' ? 'depth32float' : 'depth24plus';
    case 'depth-stencil':
      return type === 'float-stencil' ? 'depth32float-stencil8' : 'depth24plus-stencil8';
  }
  const channels = format;
  switch (type) {
    case 'ubyte':
      return `${channels}8unorm` as GPUTextureFormat;
    case 'float':
      return `${channels}32float` as GPUTextureFormat;
    case 'fp16':
      return `${channels}16float` as GPUTextureFormat;
    case 'int':
      return `${channels}32sint` as GPUTextureFormat;
    default:
      throw new Error(`unsupported texture type '${type}' for format '${format}'`);
  }
}

function getChannelCount(format: TextureFormat) {
  switch (format) {
    case 'r':
      return 1;
    case 'rg':
      return 2;
    case 'rgba':
      return 4;
    default:
      throw new Error(`no channels for format '${format}'`);
  }
}

/** Formats a filtering sampler can read; 32-bit float formats need the `float32-filterable` feature */
function isFilterableFormat(device: GPUDevice, format: GPUTextureFormat) {
  if (format.includes('depth') || format.includes('stencil') || format.endsWith('int')) return false;
  if (format.endsWith('32float')) return device.features.has('float32-filterable');
  return true;
}

/** Whether programs read `texture` with linear filtering (a linear texture of a filterable format) */
export function isLinearFiltered(device: GPUDevice, texture: WebGPUTexture) {
  return texture.filter === 'linear' && isFilterableFormat(device, texture.format);
}

export function createWebGPUTexture(
  env: WebGPUResourceEnv,
  kind: TextureKind,
  format: TextureFormat,
  type: TextureType,
  filter: TextureFilter,
): WebGPUTexture {
  const { release } = env;
  if (kind === 'texture') notSupported(`texture kind '${kind}'`);

  const dimension: GPUTextureDimension = kind.startsWith('volume-') ? '3d' : '2d';
  const gpuFormat = getGPUTextureFormat(format, type);
  const isDepth = format === 'depth' || format === 'depth-stencil';
  const hasStencil = format === 'depth-stencil';
  // all mapped color formats are renderable, 3d textures are not attachable as a whole
  const usage =
    TextureUsage.TEXTURE_BINDING |
    TextureUsage.COPY_DST |
    TextureUsage.COPY_SRC |
    (dimension === '2d' ? TextureUsage.RENDER_ATTACHMENT : 0);
  const bytesPerChannel = type === 'ubyte' ? 1 : type === 'fp16' ? 2 : 4;

  let texture: GPUTexture | undefined;
  let view: GPUTextureView | undefined;
  let attachmentView: GPUTextureView | undefined;
  let width = 0;
  let height = 0;
  let depth = 0;
  let version = 0;
  /** the data of the last `load` and whether it was mipmapped, uploaded again by `reset` */
  let loaded: Parameters<WebGPUTexture['load']>[0] | undefined;
  let mipmapped = false;

  function allocate(w: number, h: number, d = 1) {
    if (texture && w === width && h === height && d === depth) return;
    if (texture) release(texture);
    width = w;
    height = h;
    depth = d;
    texture = env.device.createTexture({
      size: { width: Math.max(1, w), height: Math.max(1, h), depthOrArrayLayers: Math.max(1, d) },
      format: gpuFormat,
      dimension,
      usage,
    });
    attachmentView = texture.createView();
    view = hasStencil ? texture.createView({ aspect: 'depth-only' }) : attachmentView;
    version += 1;
  }

  /** without data, the contents are undefined until rendered into */
  function define(w: number, h: number, d = 1) {
    loaded = undefined;
    mipmapped = false;
    allocate(w, h, d);
  }

  function load(image: TextureImage<any> | TextureVolume<any> | HTMLImageElement | ImageBitmap) {
    loaded = image;
    mipmapped = false;
    const { device } = env;
    const isElement = typeof HTMLImageElement !== 'undefined' && image instanceof HTMLImageElement;
    const isBitmap = typeof ImageBitmap !== 'undefined' && image instanceof ImageBitmap;
    if (isElement || isBitmap) {
      if (format !== 'rgba' || type !== 'ubyte') notSupported(`images in '${format}' '${type}' textures`);
      const w = isElement ? image.naturalWidth || image.width : image.width;
      const h = isElement ? image.naturalHeight || image.height : image.height;
      allocate(w, h);
      // rows top to bottom, as `texImage2D` uploads images without `UNPACK_FLIP_Y_WEBGL`
      device.queue.copyExternalImageToTexture({ source: image }, { texture: texture! }, { width: w, height: h });
      return;
    }
    if (isDepth) notSupported('loading depth textures');
    const { array, width: w, height: h } = image as TextureImage<any>;
    const d = 'depth' in image ? (image as TextureVolume<any>).depth : 1;
    if (type === 'fp16') notSupported('loading fp16 textures');
    allocate(w, h, d);

    const channels = getChannelCount(format);
    const texelSize = channels * bytesPerChannel;
    device.queue.writeTexture(
      { texture: texture! },
      (image as TextureImage<any>).flipY ? flipRows(array, w * channels, h, d) : array,
      { bytesPerRow: w * texelSize, rowsPerImage: h },
      { width: w, height: h, depthOrArrayLayers: d },
    );
  }

  const t: WebGPUTexture = {
    id: getNextTextureId(),
    format: gpuFormat,
    dimension,
    getView: () => view,
    getAttachmentView: () => attachmentView,
    getGPUTexture: () => texture,
    hasStencil,
    filter,
    get version() {
      return version;
    },

    getWidth: () => width,
    getHeight: () => height,
    getDepth: () => depth,
    getByteCount: () => width * height * depth * (channelsOf(format) * bytesPerChannel),

    define,
    load,
    mipmap: () => {
      if (!texture || dimension !== '2d' || isDepth || !isFilterableFormat(env.device, gpuFormat)) return;
      const mipLevelCount = getMipLevelCount(width, height);
      if (mipLevelCount < 2) return;
      mipmapped = true;
      // recreated with a full mip chain, level 0 copied over
      const source = texture;
      texture = env.device.createTexture({
        size: { width, height, depthOrArrayLayers: depth },
        format: gpuFormat,
        dimension,
        usage,
        mipLevelCount,
      });
      const target = texture;
      env.encode((encoder) => {
        encoder.copyTextureToTexture(
          { texture: source },
          { texture: target },
          { width, height, depthOrArrayLayers: depth },
        );
        env.generateMipmaps(encoder, target);
      });
      release(source);
      view = texture.createView();
      attachmentView = texture.createView({ baseMipLevel: 0, mipLevelCount: 1 });
      version += 1;
    },
    attachFramebuffer: (framebuffer, attachment, layer) => {
      if (layer !== undefined) notSupported('attaching texture layers');
      asWebGPUFramebuffer(framebuffer).attach(attachment, t);
    },
    detachFramebuffer: (framebuffer, attachment) => {
      asWebGPUFramebuffer(framebuffer).detach(attachment);
    },
    copyFromFramebuffer: (framebuffer, x, y, w, h, dstX, dstY) => {
      const src = asWebGPUFramebuffer(framebuffer).getColorAttachment(0);
      const srcTexture = src?.getGPUTexture();
      if (!srcTexture || !texture) throw new Error('copyFromFramebuffer needs a defined source and destination');
      // rects have a bottom-left origin (GL convention of the neutral API), textures store the top row first
      env.encode((encoder) => {
        encoder.copyTextureToTexture(
          { texture: srcTexture, origin: { x, y: src!.getHeight() - y - h } },
          { texture: texture!, origin: { x: dstX, y: height - dstY - h } },
          { width: w, height: h },
        );
      });
    },

    reset: () => {
      if (texture) release(texture);
      texture = undefined;
      view = undefined;
      attachmentView = undefined;
      const wasMipmapped = mipmapped;
      if (loaded) load(loaded);
      else if (width > 0) allocate(width, height, depth);
      if (wasMipmapped) t.mipmap();
    },
    destroy: () => {
      if (texture) release(texture);
      texture = undefined;
      view = undefined;
      attachmentView = undefined;
      loaded = undefined;
      env.restorables.delete(t);
    },
  };
  WebGPUTextures.add(t);
  env.restorables.add(t);
  return t;
}

/** Copy of `array` with the rows of each of the `depth` images in reverse order, as `UNPACK_FLIP_Y_WEBGL` uploads */
function flipRows(array: ArrayBufferView & ArrayLike<number>, rowLength: number, height: number, depth: number) {
  const source = new Uint8Array(array.buffer, array.byteOffset, array.byteLength);
  const flipped = new Uint8Array(array.byteLength);
  const rowBytes = rowLength * (array.byteLength / array.length);
  const imageBytes = rowBytes * height;
  for (let z = 0; z < depth; ++z) {
    for (let y = 0; y < height; ++y) {
      const src = z * imageBytes + y * rowBytes;
      flipped.set(source.subarray(src, src + rowBytes), z * imageBytes + (height - 1 - y) * rowBytes);
    }
  }
  return flipped;
}

function channelsOf(format: TextureFormat) {
  switch (format) {
    case 'r':
    case 'depth':
      return 1;
    case 'rg':
    case 'depth-stencil':
      return 2;
    default:
      return 4;
  }
}

/** Cube faces in the array layer order of cube views (+x, -x, +y, -y, +z, -z) */
export const CubeLayerSides = ['px', 'nx', 'py', 'ny', 'pz', 'nz'] as const;

/**
 * Cube texture (rgba, linear) of six square faces, sampled through a cube view; `loadFaces` takes the face images in
 * `CubeLayerSides` order, rows top to bottom as on WebGL
 */
export function createWebGPUCubeTexture(env: WebGPUResourceEnv, mipmaps: boolean) {
  const { release } = env;
  const gpuFormat: GPUTextureFormat = 'rgba8unorm';
  let texture: GPUTexture | undefined;
  let view: GPUTextureView | undefined;
  let size = 0;
  let version = 0;

  /** the faces, loaded again by `reset` */
  let faces: ImageBitmap[] | undefined;

  function loadFaces(images: ImageBitmap[]) {
    const { device } = env;
    faces = images;
    if (texture) release(texture);
    size = images[0].width;
    const mipLevelCount = mipmaps ? getMipLevelCount(size, size) : 1;
    texture = device.createTexture({
      size: { width: size, height: size, depthOrArrayLayers: 6 },
      format: gpuFormat,
      usage: TextureUsage.TEXTURE_BINDING | TextureUsage.COPY_DST | TextureUsage.RENDER_ATTACHMENT,
      mipLevelCount,
      label: 'cube',
    });
    for (let i = 0; i < 6; ++i) {
      device.queue.copyExternalImageToTexture(
        { source: images[i] },
        { texture, origin: { x: 0, y: 0, z: i } },
        { width: size, height: size },
      );
    }
    const target = texture;
    if (mipLevelCount > 1) env.encode((encoder) => env.generateMipmaps(encoder, target));
    view = texture.createView({ dimension: 'cube' });
    version += 1;
  }

  const t: WebGPUTexture & { loadFaces: typeof loadFaces } = {
    id: getNextTextureId(),
    format: gpuFormat,
    dimension: '2d',
    getView: () => view,
    getAttachmentView: () => undefined,
    getGPUTexture: () => texture,
    hasStencil: false,
    filter: 'linear',
    get version() {
      return version;
    },

    getWidth: () => size,
    getHeight: () => size,
    getDepth: () => 0,
    getByteCount: () => size * size * 4 * 6,

    define: () => notSupported('defining cube textures'),
    load: () => notSupported('loading cube textures, use loadFaces'),
    loadFaces,
    mipmap: () => {},
    attachFramebuffer: () => notSupported('attaching cube textures'),
    detachFramebuffer: () => notSupported('detaching cube textures'),
    copyFromFramebuffer: () => notSupported('copying into cube textures'),

    reset: () => {
      if (faces) loadFaces(faces);
    },
    destroy: () => {
      if (texture) release(texture);
      texture = undefined;
      view = undefined;
      faces = undefined;
      env.restorables.delete(t);
    },
  };
  WebGPUTextures.add(t);
  env.restorables.add(t);
  return t;
}
