/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 *
 * ported from https://github.com/photopea/UZIP.js/blob/master/UZIP.js
 * MIT License, Copyright (c) 2018 Photopea
 */

import { InflateContext, inflateBlocks } from '../internal/inflate-blocks.js';
import { checkCompressionStreamSupport } from './util.js';
import type { RuntimeContext } from '@molstar/core/task';

// inflating a 44 MB gzip file to 225 MB.

// Using JS
// inflate: 570.8759765625 ms
// inflate: 562.950927734375 ms
// inflate: 595.06591796875 ms
// inflate: 584.8740234375 ms
// inflate: 583.31396484375 ms

// Using DecompressionStream
// inflate: 502.001953125 ms
// inflate: 476.818115234375 ms
// inflate: 476.68701171875 ms
// inflate: 422.319091796875 ms
// inflate: 428.925048828125 ms

/**
 * Inflates `data` with the platform's `DecompressionStream`, into `buf` when given (never past it: a longer stream
 * throws), and returns what it inflated. `'gzip'` checks the member's CRC-32 and size itself.
 */
export async function inflateStream(
  runtime: RuntimeContext,
  data: Uint8Array<ArrayBuffer>,
  format: 'deflate-raw' | 'gzip',
  buf?: Uint8Array<ArrayBuffer>,
): Promise<Uint8Array<ArrayBuffer>> {
  const ds = new DecompressionStream(format);
  const blob = new Blob([data]);
  const decompressedStream = blob.stream().pipeThrough(ds);

  let offset = 0;
  const chunks: Uint8Array[] = [];
  const reader = decompressedStream.getReader();

  const readChunk = async (): Promise<undefined> => {
    const { done, value } = await reader.read();
    if (done) return;
    // Never past the buffer given (the declared size): stopped, not grown
    if (buf && offset + value.length > buf.length) {
      await reader.cancel();
      throw new Error(`inflates past its declared ${buf.length} bytes`);
    }

    if (runtime.shouldUpdate) {
      await runtime.update({ message: 'Inflating blocks...', current: offset, max: buf?.length });
    }
    if (buf) {
      buf.set(value, offset);
    } else {
      chunks.push(value);
    }
    offset += value.length;
    return readChunk();
  };
  await readChunk();

  if (!buf) {
    buf = new Uint8Array(offset);
    for (let i = 0, j = 0; i < chunks.length; i++) {
      buf.set(chunks[i], j);
      j += chunks[i].length;
    }
  }
  // What was inflated, so a caller sees a stream shorter than the buffer
  return offset < buf.length ? buf.subarray(0, offset) : buf;
}

// https://tools.ietf.org/html/rfc1951
export async function _inflate(
  runtime: RuntimeContext,
  data: Uint8Array<ArrayBuffer>,
  buf?: Uint8Array<ArrayBuffer>,
): Promise<Uint8Array<ArrayBuffer>> {
  if (data[0] === 3 && data[1] === 0) return buf ? buf : new Uint8Array(0);

  if (checkCompressionStreamSupport('deflate-raw')) return inflateStream(runtime, data, 'deflate-raw', buf);

  const ctx = InflateContext(data, buf);
  while (ctx.BFINAL === 0) {
    if (runtime.shouldUpdate) {
      await runtime.update({ message: 'Inflating blocks...', current: ctx.pos, max: data.length });
    }
    inflateBlocks(ctx, 100);
  }
  return ctx.buf.length === ctx.off ? ctx.buf : ctx.buf.slice(0, ctx.off);
}
