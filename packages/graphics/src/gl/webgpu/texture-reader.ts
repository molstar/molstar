/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { isDebugMode } from '@molstar/core/util/debug';
import type { Texture } from '../texture.js';
import type {
  TextureReadArray,
  TextureReader,
  TextureReaderOptions,
  TextureReadFormat,
  TextureReadStatus,
} from '../texture-reader.js';
import { BufferUsage, MapMode } from './constants.js';

/** What a reader needs of its context */
export interface WebGPUTextureReadContext {
  readonly device: GPUDevice;
  /** submits the commands recorded so far: a read is ordered after them */
  submit: () => void;
  /** the texture of `texture` to read as `format`; throws where it has no storage or another format */
  getTexture: (texture: Texture, format: TextureReadFormat) => GPUTexture;
  /** the staging buffers' bytes and count added to the context's (negative: taken) */
  countStaging: (bytes: number, buffers: number) => void;
}

/**
 * Reads asynchronously only: a read's textures copied into one staging buffer in one command buffer, mapped once, its
 * rows copied straight into the caller's arrays; the buffer kept for the next read (a second one while a read waits)
 */
export function createWebGPUTextureReader<F extends TextureReadFormat = 'rgba8'>(
  ctx: WebGPUTextureReadContext,
  options: TextureReaderOptions<F> = {},
): TextureReader<TextureReadArray<F>> {
  const format = options.format ?? 'rgba8';
  let status: TextureReadStatus = 'failed';
  /** results of reads started before the last one (or a reset) are dropped */
  let readId = 0;
  /** the buffer kept for the next read, and the one of the read resolved (mapped) */
  let spare: GPUBuffer | undefined;
  let mapped: GPUBuffer | undefined;
  let layout = { count: 0, rowSize: 0, bytesPerRow: 0, height: 0 };
  let disposed = false;

  function release(buffer: GPUBuffer) {
    ctx.countStaging(-buffer.size, -1);
    buffer.destroy();
  }

  /** `buffer` kept as the spare (the larger of two), unmapped */
  function keep(buffer: GPUBuffer) {
    if (buffer.mapState === 'mapped') buffer.unmap();
    if (disposed || (spare && spare.size >= buffer.size)) release(buffer);
    else {
      if (spare) release(spare);
      spare = buffer;
    }
  }

  function drop() {
    readId += 1;
    status = 'failed';
    if (mapped) keep(mapped);
    mapped = undefined;
  }

  return {
    synchronous: false,
    asynchronous: true,
    read: () => {
      throw new Error('synchronous texture reads are not supported');
    },
    start: (textures, x, y, width, height) => {
      const gpuTextures = textures.map((t) => ctx.getTexture(t, format));
      drop();
      const id = readId;
      status = 'pending';
      ctx.submit();
      // four bytes a texel (`TextureReadFormat`), rows at 256 byte strides
      const rowSize = width * 4;
      const bytesPerRow = Math.ceil(rowSize / 256) * 256;
      const size = bytesPerRow * height * gpuTextures.length;
      let buffer = spare;
      spare = undefined;
      if (!buffer || buffer.size < size) {
        if (buffer) release(buffer);
        buffer = ctx.device.createBuffer({ size, usage: BufferUsage.COPY_DST | BufferUsage.MAP_READ });
        ctx.countStaging(size, 1);
      }
      const encoder = ctx.device.createCommandEncoder({ label: 'texture read' });
      for (let i = 0; i < gpuTextures.length; ++i) {
        const offset = i * bytesPerRow * height;
        encoder.copyTextureToBuffer(
          { texture: gpuTextures[i], origin: { x, y } },
          { buffer, offset, bytesPerRow },
          { width, height },
        );
      }
      ctx.device.queue.submit([encoder.finish()]);
      const b = buffer;
      b.mapAsync(MapMode.READ, 0, size).then(
        () => {
          if (id !== readId) return keep(b);
          mapped = b;
          layout = { count: gpuTextures.length, rowSize, bytesPerRow, height };
          status = 'resolved';
        },
        (e) => {
          // the reader disposed (its buffer destroyed) or the device lost
          if (id !== readId) return;
          if (isDebugMode) console.error(e);
          status = 'failed';
        },
      );
    },
    check: (arrays) => {
      if (status !== 'resolved' || !mapped) return status;
      const { count, rowSize, bytesPerRow, height } = layout;
      const src = new Uint8Array(mapped.getMappedRange(0, bytesPerRow * height * count));
      // rows top to bottom, as read; each array takes what fits
      for (let i = 0; i < Math.min(arrays.length, count); ++i) {
        const target = new Uint8Array(arrays[i].buffer, arrays[i].byteOffset, arrays[i].byteLength);
        for (let r = 0, o = 0; r < height && o < target.length; ++r, o += rowSize) {
          const start = (i * height + r) * bytesPerRow;
          target.set(src.subarray(start, start + Math.min(rowSize, target.length - o)), o);
        }
      }
      drop();
      return 'resolved';
    },
    reset: drop,
    getByteCount: () => (spare?.size ?? 0) + (mapped?.size ?? 0),
    dispose: () => {
      drop();
      disposed = true;
      if (spare) release(spare);
      spare = undefined;
    },
  };
}
