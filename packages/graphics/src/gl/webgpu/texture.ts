/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { toHalfFloat } from '@molstar/core/util/number-conversion';
import { idFactory } from '@molstar/core/util/id-factory';
import type { TextureImage, TextureVolume } from '../renderable/util.js';
import type { GraphicsStats } from '../stats.js';
import type { Texture, TextureFilter, TextureFormat, TextureKind, TextureSource, TextureType } from '../texture.js';
import { TextureUsage } from './constants.js';
import { notSupported } from './errors.js';
import type { AttachmentRef } from '../framebuffer.js';
import { asWebGPUFramebuffer } from './framebuffer.js';
import { getMipLevelCount, type MipmapGenerator } from './mipmap.js';

const getNextTextureId = idFactory();

/**
 * Destroys a GPU resource once no recorded commands use it: WebGPU invalidates destroyed resources for
 * all later submits, also for commands recorded before
 */
export type Release = (resource: { destroy(): void }) => void;

/** What WebGPU resources need from their context */
/** Bytes of the live resources of a context by kind, as `GraphicsResources.getByteCounts` reports them */
export type WebGPUByteCounts = {
  texture: number;
  cubeTexture: number;
  attribute: number;
  elements: number;
  pixelPack: number;
  uniform: number;
  renderbuffer: number;
};

export function createWebGPUByteCounts(): WebGPUByteCounts {
  return { texture: 0, cubeTexture: 0, attribute: 0, elements: 0, pixelPack: 0, uniform: 0, renderbuffer: 0 };
}

/**
 * Keeps the bytes of one resource in `counts[kind]`: `set` replaces them (on allocation, 0 when released). Each
 * resource accounts its own bytes, so that a total is the sum of what is live, as WebGL sums its resources.
 */
export function createByteAccount(counts: WebGPUByteCounts, kind: keyof WebGPUByteCounts) {
  let bytes = 0;
  return (next: number) => {
    counts[kind] += next - bytes;
    bytes = next;
  };
}

/**
 * Counts a resource in `stats.resourceCounts[name]` from now until the returned function is called (once, on destroy),
 * as WebGL counts its live resources
 */
export function countResource(stats: GraphicsStats, name: keyof GraphicsStats['resourceCounts']) {
  stats.resourceCounts[name] += 1;
  let counted = true;
  return () => {
    if (!counted) return;
    counted = false;
    stats.resourceCounts[name] -= 1;
  };
}

export interface WebGPUResourceEnv {
  readonly device: GPUDevice;
  /** counts of the live resources (`resourceCounts`) */
  readonly stats: GraphicsStats;
  /** bytes of the live resources, see `WebGPUByteCounts` */
  readonly byteCounts: WebGPUByteCounts;
  readonly release: Release;
  /** Encodes commands outside of render passes (copies), ends the open pass */
  readonly encode: (commands: (encoder: GPUCommandEncoder) => void) => void;
  /** Submits the recorded commands */
  readonly flush: () => void;
  readonly generateMipmaps: MipmapGenerator;
}

export interface WebGPUTexture extends Texture {
  /** as attachment (`Framebuffer.getAttachments`) */
  readonly ref: AttachmentRef;
  readonly format: GPUTextureFormat;
  readonly dimension: GPUTextureDimension;
  /**
   * view for sampling (depth aspect of depth-stencil formats), undefined while the texture has no storage: before it is
   * first written (bindings use a placeholder then, see `getTextureEntries`)
   */
  getView: () => GPUTextureView | undefined;
  /** view for render pass attachments and storage bindings (all aspects): allocates the storage of a defined texture */
  getAttachmentView: () => GPUTextureView | undefined;
  /** the GPU texture, for copies and reads: allocates the storage of a defined texture */
  getGPUTexture: () => GPUTexture | undefined;
  /** whether the texture has GPU storage (written, or read or copied since `define`) */
  readonly hasStorage: boolean;
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

/** `features` of the device: `depth32float-stencil8` is optional, `depth24plus-stencil8` is used without it */
export function getGPUTextureFormat(
  format: TextureFormat,
  type: TextureType,
  features?: GPUSupportedFeatures,
): GPUTextureFormat {
  switch (format) {
    case 'depth':
      return type === 'float' ? 'depth32float' : 'depth24plus';
    case 'depth-stencil':
      return type === 'float-stencil' && features?.has('depth32float-stencil8')
        ? 'depth32float-stencil8'
        : 'depth24plus-stencil8';
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

/** Formats with storage texture support in core WebGPU */
const StorageFormats = new Set<GPUTextureFormat>([
  'rgba8unorm',
  'rgba16float',
  'rgba32float',
  'r32float',
  'r32sint',
  'r32uint',
  'rg32float',
]);

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
  const gpuFormat = getGPUTextureFormat(format, type, env.device.features);
  const isDepth = format === 'depth' || format === 'depth-stencil';
  const hasStencil = format === 'depth-stencil';
  // all mapped color formats are renderable, 3d textures one slice at a time (`attachFramebuffer` with a layer); formats
  // with storage support can be written by compute kernels
  const usage =
    TextureUsage.TEXTURE_BINDING |
    TextureUsage.COPY_DST |
    TextureUsage.COPY_SRC |
    TextureUsage.RENDER_ATTACHMENT |
    (StorageFormats.has(gpuFormat) ? TextureUsage.STORAGE_BINDING : 0);
  const bytesPerChannel = type === 'ubyte' ? 1 : type === 'fp16' ? 2 : 4;
  const account = createByteAccount(env.byteCounts, 'texture');
  const uncount = countResource(env.stats, 'texture');

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
  /** a destroyed texture allocates nothing again (e.g. an image decoded after the background changed) */
  let destroyed = false;
  /** defined (`define`) without storage yet: allocated by the first write, copy or read (`ensureStorage`) */
  let defined = false;

  function allocate(w: number, h: number, d = 1) {
    if (texture && w === width && h === height && d === depth) return;
    if (texture) release(texture);
    width = w;
    height = h;
    depth = d;
    defined = false;
    account(Math.max(1, w) * Math.max(1, h) * Math.max(1, d) * channelsOf(format) * bytesPerChannel);
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

  /** `allocate` for a queue write */
  function allocateForWrite(w: number, h: number, d = 1) {
    // queue writes run before the commands recorded so far, which may still read the current contents (uploads are
    // ordered on WebGL): submit them first
    if (texture && w === width && h === height && d === depth) env.flush();
    allocate(w, h, d);
  }

  /**
   * Records the size; the storage is allocated when the texture is first written, copied or read (lazy: targets of
   * passes that never draw take no memory). Without data, the contents are undefined until rendered into.
   */
  function define(w: number, h: number, d = 1) {
    if (destroyed) return;
    loaded = undefined;
    mipmapped = false;
    // storage of this size is kept (its contents are undefined either way)
    if (texture && w === width && h === height && d === depth) return;
    if (texture) {
      release(texture);
      texture = undefined;
      view = undefined;
      attachmentView = undefined;
      account(0);
      version += 1;
    }
    width = w;
    height = h;
    depth = d;
    defined = true;
  }

  /** Allocates the storage of a defined texture, see `define` */
  function ensureStorage() {
    if (!texture && defined && !destroyed) allocate(width, height, depth);
  }

  function load(image: TextureImage<any> | TextureVolume<any> | HTMLImageElement | ImageBitmap) {
    if (destroyed) return;
    loaded = image;
    mipmapped = false;
    const { device } = env;
    const isElement = typeof HTMLImageElement !== 'undefined' && image instanceof HTMLImageElement;
    const isBitmap = typeof ImageBitmap !== 'undefined' && image instanceof ImageBitmap;
    if (isElement || isBitmap) {
      if (format !== 'rgba' || type !== 'ubyte') notSupported(`images in '${format}' '${type}' textures`);
      const w = isElement ? image.naturalWidth || image.width : image.width;
      const h = isElement ? image.naturalHeight || image.height : image.height;
      allocateForWrite(w, h);
      // rows top to bottom, as `texImage2D` uploads images without `UNPACK_FLIP_Y_WEBGL`
      device.queue.copyExternalImageToTexture({ source: image }, { texture: texture! }, { width: w, height: h });
      return;
    }
    if (isDepth) notSupported('loading depth textures');
    const { width: w, height: h } = image as TextureImage<any>;
    let { array } = image as TextureImage<any>;
    const d = 'depth' in image ? (image as TextureVolume<any>).depth : 1;
    // half floats as their bits (`Uint16Array`, as WebGL's `HALF_FLOAT` takes them); floats are converted
    if (type === 'fp16' && array instanceof Float32Array) array = Uint16Array.from(array, toHalfFloat);
    allocateForWrite(w, h, d);

    const channels = getChannelCount(format);
    const texelSize = channels * bytesPerChannel;
    device.queue.writeTexture(
      { texture: texture! },
      (image as TextureImage<any>).flipY ? flipRows(array, w * channels, h, d) : array,
      { bytesPerRow: w * texelSize, rowsPerImage: h },
      { width: w, height: h, depthOrArrayLayers: d },
    );
  }

  /** `Texture.loadRows`: a queue write per range of rows, from the array in place */
  function loadRows(image: TextureImage<any>, rows: readonly number[]) {
    if (destroyed) return;
    const { array } = image;
    const converted = type === 'fp16' && array instanceof Float32Array;
    if (
      dimension !== '2d' ||
      image.flipY ||
      converted ||
      !texture ||
      image.width !== width ||
      image.height !== height
    ) {
      load(image);
      return;
    }
    // queue writes run before the commands recorded so far, which may still read the current contents
    env.flush();
    const stride = (array.length / image.height) * array.BYTES_PER_ELEMENT;
    for (let k = 0; k < rows.length; k += 2) {
      const start = rows[k];
      const end = Math.min(rows[k + 1], image.height);
      if (end <= start) continue;
      env.device.queue.writeTexture(
        { texture, origin: { x: 0, y: start } },
        array,
        { offset: start * stride, bytesPerRow: stride, rowsPerImage: end - start },
        { width, height: end - start },
      );
    }
    loaded = image;
  }

  const id = getNextTextureId();
  const t: WebGPUTexture = {
    id,
    ref: { kind: 'texture', id },
    format: gpuFormat,
    dimension,
    getView: () => view,
    getAttachmentView: () => {
      ensureStorage();
      return attachmentView;
    },
    getGPUTexture: () => {
      ensureStorage();
      return texture;
    },
    get hasStorage() {
      return texture !== undefined;
    },
    getSource: () => ({ kind, format, type, filter, data: loaded as TextureSource['data'], mipmap: mipmapped }),
    hasStencil,
    filter,
    get version() {
      return version;
    },

    getWidth: () => width,
    getHeight: () => height,
    getDepth: () => depth,
    getByteCount: () => (texture ? width * height * depth * (channelsOf(format) * bytesPerChannel) : 0),

    define,
    load,
    loadRows,
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
      if (layer === undefined) {
        asWebGPUFramebuffer(framebuffer).attach(attachment, t);
        return;
      }
      if (dimension !== '3d') notSupported('attaching layers of 2d textures');
      // a slice of the 3d view of level 0
      asWebGPUFramebuffer(framebuffer).attach(attachment, {
        ref: t.ref,
        format: gpuFormat,
        hasStencil,
        depthSlice: layer,
        getAttachmentView: () => {
          ensureStorage();
          return attachmentView;
        },
        getGPUTexture: () => {
          ensureStorage();
          return texture;
        },
        getWidth: () => width,
        getHeight: () => height,
      });
    },
    detachFramebuffer: (framebuffer, attachment) => {
      asWebGPUFramebuffer(framebuffer).detach(attachment);
    },
    copyFromFramebuffer: (framebuffer, x, y, w, h, dstX, dstY) => {
      const src = asWebGPUFramebuffer(framebuffer).getColorAttachment(0);
      const srcTexture = src?.getGPUTexture();
      ensureStorage();
      if (!srcTexture || !texture) throw new Error('copyFromFramebuffer needs a defined source and destination');
      env.encode((encoder) => {
        encoder.copyTextureToTexture(
          { texture: srcTexture, origin: { x, y } },
          { texture: texture!, origin: { x: dstX, y: dstY } },
          { width: w, height: h },
        );
      });
    },

    releaseStorage: () => {
      if (!texture || destroyed) return;
      release(texture);
      texture = undefined;
      view = undefined;
      attachmentView = undefined;
      loaded = undefined;
      mipmapped = false;
      account(0);
      version += 1;
      // allocated again at this size by the next write
      defined = true;
    },
    destroy: () => {
      if (destroyed) return;
      destroyed = true;
      if (texture) release(texture);
      texture = undefined;
      view = undefined;
      attachmentView = undefined;
      loaded = undefined;
      defined = false;
      account(0);
      uncount();
    },
  };
  WebGPUTextures.add(t);
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
  const account = createByteAccount(env.byteCounts, 'cubeTexture');
  const uncount = countResource(env.stats, 'cubeTexture');
  let texture: GPUTexture | undefined;
  let view: GPUTextureView | undefined;
  let size = 0;
  let version = 0;
  /** a destroyed texture allocates nothing again (faces decoded after the background changed) */
  let destroyed = false;

  function loadFaces(images: ImageBitmap[]) {
    if (destroyed) return;
    const { device } = env;
    if (texture) release(texture);
    size = images[0].width;
    account(size * size * 4 * 6);
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

  const cubeId = getNextTextureId();
  const t: WebGPUTexture & { loadFaces: typeof loadFaces } = {
    id: cubeId,
    ref: { kind: 'texture', id: cubeId },
    format: gpuFormat,
    dimension: '2d',
    getView: () => view,
    getAttachmentView: () => undefined,
    getGPUTexture: () => texture,
    get hasStorage() {
      return texture !== undefined;
    },
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
    // loaded from images, not a target
    releaseStorage: () => {},

    destroy: () => {
      if (destroyed) return;
      destroyed = true;
      if (texture) release(texture);
      texture = undefined;
      view = undefined;
      account(0);
      uncount();
    },
  };
  WebGPUTextures.add(t);
  return t;
}
