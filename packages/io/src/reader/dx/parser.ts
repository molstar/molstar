/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from NGL.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { chunkedSubtask, type RuntimeContext, Task } from '@molstar/core/task';
import type { StringLike } from '@molstar/core/util/string-like';
import { parseFloatBytes } from '@molstar/core/util/number-parser';
import { TextByteTokenizer } from '../common/text/text-byte-tokenizer.js';
import { ReaderResult as Result } from '../result.js';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';

// http://apbs-pdb2pqr.readthedocs.io/en/latest/formats/opendx.html

export interface DxFile {
  name: string;
  header: DxFile.Header;
  values: Float64Array;
}

export namespace DxFile {
  export interface Header {
    dim: Vec3;
    min: Vec3;
    h: Vec3;
  }
}

/** Reads the header lines (decoded, few) up to `object 3`; the tokenizer is then at the first byte of the values. */
function readHeader(tokenizer: TextByteTokenizer): DxFile.Header {
  const header: Partial<DxFile.Header> = { h: Vec3() };
  let deltaLineCount = 0;

  const reWhitespace = /\s+/g;

  while (true) {
    if (tokenizer.position >= tokenizer.length) throw new Error(`no 'object 3' record found`);
    const line = TextByteTokenizer.readLine(tokenizer);
    let ls;

    if (line.startsWith('object 1')) {
      ls = line.split(reWhitespace);
      header.dim = Vec3.create(parseInt(ls[5]), parseInt(ls[6]), parseInt(ls[7]));
    } else if (line.startsWith('origin')) {
      ls = line.split(reWhitespace);
      header.min = Vec3.create(parseFloat(ls[1]), parseFloat(ls[2]), parseFloat(ls[3]));
    } else if (line.startsWith('delta')) {
      ls = line.split(reWhitespace);

      if (deltaLineCount === 0) {
        (header.h as any)[0] = parseFloat(ls[1]);
      } else if (deltaLineCount === 1) {
        (header.h as any)[1] = parseFloat(ls[2]);
      } else if (deltaLineCount === 2) {
        (header.h as any)[2] = parseFloat(ls[3]);
      }

      deltaLineCount += 1;
    } else if (line.startsWith('object 3')) {
      break;
    }
  }

  return header as DxFile.Header;
}

function readValuesText(ctx: RuntimeContext, tokenizer: TextByteTokenizer, header: DxFile.Header) {
  const N = header.dim[0] * header.dim[1] * header.dim[2];
  const chunkSize = 100 * 100 * 100;
  const values = new Float64Array(N);
  const { data, length } = tokenizer;
  let offset = 0;

  return chunkedSubtask(
    ctx,
    chunkSize,
    values,
    (count, values) => {
      const max = Math.min(N, offset + count);
      // locals, not the tokenizer's fields, per byte: the hot loop of a large grid
      let p = tokenizer.position;
      for (let i = offset; i < max; i++) {
        let c = data[p];
        while (p < length && (c === 32 || c === 10 || c === 13 || c === 9)) c = data[++p];
        const start = p;
        while (p < length && c !== 32 && c !== 10 && c !== 13 && c !== 9) c = data[++p];
        values[i] = parseFloatBytes(data, start, p);
      }
      tokenizer.position = p;
      offset = max;
      return max === N ? 0 : chunkSize;
    },
    (ctx, _, i) => ctx.update({ current: Math.min(i, N), max: N }),
  );
}

async function parseText(taskCtx: RuntimeContext, data: ReadonlyUint8Array, name: string) {
  await taskCtx.update('Reading header...');
  const tokenizer = TextByteTokenizer(data);
  const header = readHeader(tokenizer);
  await taskCtx.update('Reading values...');
  const values = await readValuesText(taskCtx, tokenizer, header);
  return Result.success({ header, values, name });
}

async function parseBinary(taskCtx: RuntimeContext, data: ReadonlyUint8Array, name: string) {
  await taskCtx.update('Reading header...');

  const tokenizer = TextByteTokenizer(data);
  const header = readHeader(tokenizer);
  // the values start right after the `object 3` line
  const headerByteCount = tokenizer.position;

  await taskCtx.update('Reading values...');

  const size = header.dim[0] * header.dim[1] * header.dim[2];
  const dv = new DataView(data.buffer, data.byteOffset + headerByteCount);
  const values = new Float64Array(size);

  for (let i = 0; i < size; i++) {
    values[i] = dv.getFloat64(i * 8, true);
  }

  // TODO: why doesnt this work? throw "attempting to construct out-of-bounds TypedArray"
  // const values = new Float64Array(data.buffer, data.byteOffset + headerByteCount, header.dim[0] * header.dim[1] * header.dim[2]);
  return Result.success({ header, values, name });
}

/**
 * Parses text DX from its text or its UTF-8 bytes, or binary DX (`.dxbin`: the text header, then little-endian
 * float64 values) from its bytes. `binary` defaults to whether `data` is bytes; a caller with the bytes of a text DX
 * passes `binary: false`.
 */
export function parseDx(data: StringLike | ReadonlyUint8Array, name: string, options?: { binary?: boolean }) {
  return Task.create<Result<DxFile>>('Parse DX', (taskCtx) => {
    const binary = options?.binary ?? data instanceof Uint8Array;
    const bytes = TextByteTokenizer.bytesOf(data);
    if (binary) return parseBinary(taskCtx, bytes, name);
    return parseText(taskCtx, bytes, name);
  });
}
