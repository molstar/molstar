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

/** Reads `width * height` texels at (x, y) of `texture` (a top-left origin), rows top to bottom */
export type WebGPUTextureRead = (
  texture: Texture,
  format: TextureReadFormat,
  x: number,
  y: number,
  width: number,
  height: number,
) => Promise<Uint8Array>;

/** Reads asynchronously only, through staging buffers mapped once the copies completed */
export function createWebGPUTextureReader<F extends TextureReadFormat = 'rgba8'>(
  readTexture: WebGPUTextureRead,
  options: TextureReaderOptions<F> = {},
): TextureReader<TextureReadArray<F>> {
  const format = options.format ?? 'rgba8';
  let status: TextureReadStatus = 'failed';
  let data: Uint8Array[] = [];
  /** results of reads started before the last one (or a reset) are dropped */
  let readId = 0;

  function drop() {
    readId += 1;
    status = 'failed';
    data = [];
  }

  return {
    synchronous: false,
    asynchronous: true,
    read: () => {
      throw new Error('synchronous texture reads are not supported');
    },
    start: (textures, x, y, width, height) => {
      const id = ++readId;
      status = 'pending';
      Promise.all(textures.map((t) => readTexture(t, format, x, y, width, height))).then(
        (arrays) => {
          if (id !== readId) return;
          data = arrays;
          status = 'resolved';
        },
        (e) => {
          if (id !== readId) return;
          if (isDebugMode) console.error(e);
          status = 'failed';
        },
      );
    },
    check: (arrays) => {
      if (status !== 'resolved') return status;
      // rows top to bottom, as read
      for (let i = 0; i < arrays.length; ++i) {
        const target = new Uint8Array(arrays[i].buffer, arrays[i].byteOffset, arrays[i].byteLength);
        target.set(data[i].subarray(0, target.length));
      }
      drop();
      return 'resolved';
    },
    reset: drop,
    // staging buffers live only while a read is pending
    getByteCount: () => 0,
    dispose: drop,
  };
}
