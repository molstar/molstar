/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Texture } from '../texture.js';
import type {
  TextureReadArray,
  TextureReader,
  TextureReaderOptions,
  TextureReadFormat,
  TextureReadStatus,
} from '../texture-reader.js';
import type { PixelPackBuffer } from './buffer.js';
import { isWebGL2 } from './compat.js';
import type { WebGLContext } from './context.js';
import { endRenderPass } from './render-pass.js';

/** Reverses the order of the `rows` rows of `rowLength` values of `array`, GL reads rows bottom to top */
function flipRows(array: Uint8Array | Float32Array, rowLength: number, rows: number) {
  const tmp = array.slice(0, rowLength);
  for (let a = 0, b = rows - 1; a < b; ++a, --b) {
    const ia = a * rowLength;
    const ib = b * rowLength;
    tmp.set(array.subarray(ia, ia + rowLength));
    array.copyWithin(ia, ib, ib + rowLength);
    array.set(tmp, ib);
  }
}

/**
 * Reads with `readPixels`; asynchronously on WebGL2 through pixel pack buffers and a fence, which `check` polls.
 * `r32f` reads red/float where the implementation allows it, else rgba/float (always readable from float color buffers)
 * keeping the red channel.
 */
export function createWebGLTextureReader<F extends TextureReadFormat = 'rgba8'>(
  ctx: WebGLContext,
  options: TextureReaderOptions<F> = {},
): TextureReader<TextureReadArray<F>> {
  const { gl, resources, state } = ctx;
  const isFloat = options.format === 'r32f';
  const framebuffer = resources.framebuffer();
  const buffers: PixelPackBuffer[] = [];

  /** for `r32f`, whether red/float is readable (else rgba/float), known once a float texture was bound */
  let readsRed: boolean | undefined;
  /** rgba/float texels of `r32f` reads without red/float reads */
  let rgbaFloats = new Float32Array(0);
  let fenceSync: WebGLSync | null = null;
  let byteCount = 0;
  /** the size of the rectangle of the started read */
  let readWidth = 0;
  let readHeight = 0;
  /** values per texel of the result arrays */
  const texelLength = isFloat ? 1 : 4;

  /** binds the read framebuffer with `texture` attached */
  function bind(texture: Texture) {
    texture.attachFramebuffer(framebuffer, 'color0');
    framebuffer.bind();
    if (isFloat && readsRed === undefined) {
      readsRed =
        gl.getParameter(gl.IMPLEMENTATION_COLOR_READ_FORMAT) === (gl as WebGL2RenderingContext).RED &&
        gl.getParameter(gl.IMPLEMENTATION_COLOR_READ_TYPE) === gl.FLOAT;
    }
  }

  /** an array for rgba/float reads of `count` texels */
  function getRgbaFloats(count: number) {
    if (rgbaFloats.length < count * 4) rgbaFloats = new Float32Array(count * 4);
    return rgbaFloats.subarray(0, count * 4);
  }

  function copyRed(rgba: Float32Array, array: Float32Array, count: number) {
    for (let i = 0; i < count; ++i) array[i] = rgba[i * 4];
  }

  function deleteFence() {
    if (fenceSync !== null && isWebGL2(gl)) gl.deleteSync(fenceSync);
    fenceSync = null;
  }

  return {
    synchronous: true,
    asynchronous: isWebGL2(gl),
    read: (textures, x, topY, width, height, arrays) => {
      endRenderPass(state);
      for (let i = 0; i < textures.length; ++i) {
        bind(textures[i]);
        // the rectangle has a top-left origin, GL's a bottom-left one
        const y = textures[i].getHeight() - topY - height;
        if (!isFloat) {
          ctx.readPixels(x, y, width, height, arrays[i]);
        } else if (readsRed && isWebGL2(gl)) {
          gl.readPixels(x, y, width, height, gl.RED, gl.FLOAT, arrays[i]);
        } else {
          const rgba = getRgbaFloats(width * height);
          ctx.readPixels(x, y, width, height, rgba);
          copyRed(rgba, arrays[i] as Float32Array, width * height);
        }
        flipRows(arrays[i], width * texelLength, height);
      }
    },
    start: (textures, x, topY, width, height) => {
      if (!isWebGL2(gl)) return;
      endRenderPass(state);
      deleteFence();
      byteCount = 0;
      readWidth = width;
      readHeight = height;
      for (let i = 0; i < textures.length; ++i) {
        bind(textures[i]);
        if (!buffers[i]) {
          buffers[i] = isFloat
            ? resources.pixelPack(readsRed ? 'r' : 'rgba', 'float')
            : resources.pixelPack('rgba', 'ubyte');
        }
        // the rectangle has a top-left origin, GL's a bottom-left one
        buffers[i].read(x, textures[i].getHeight() - topY - height, width, height);
        byteCount += buffers[i].getByteCount();
      }
      fenceSync = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
      gl.flush();
    },
    check: (arrays): TextureReadStatus => {
      if (fenceSync === null || !isWebGL2(gl)) return 'failed';
      const res = gl.clientWaitSync(fenceSync, 0, 0);
      if (res === gl.WAIT_FAILED) {
        deleteFence();
        return 'failed';
      }
      if (res === gl.TIMEOUT_EXPIRED) return 'pending';
      for (let i = 0; i < arrays.length; ++i) {
        if (isFloat && !readsRed) {
          const count = arrays[i].length;
          const rgba = getRgbaFloats(count);
          buffers[i].getSubData(rgba);
          copyRed(rgba, arrays[i] as Float32Array, count);
        } else {
          buffers[i].getSubData(arrays[i]);
        }
        flipRows(arrays[i], readWidth * texelLength, readHeight);
      }
      deleteFence();
      return 'resolved';
    },
    reset: () => {
      deleteFence();
    },
    getByteCount: () => byteCount,
    dispose: () => {
      deleteFence();
      for (const b of buffers) b.destroy();
      buffers.length = 0;
      byteCount = 0;
      framebuffer.destroy();
    },
  };
}
