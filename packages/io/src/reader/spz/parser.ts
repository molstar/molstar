/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type RuntimeContext, Task } from '@molstar/core/task';
import type { ReadonlyUint8Array } from '@molstar/core/util/type-helpers';
import { isGzip, ungzip } from '@molstar/core/util/zip/zip';
import { ReaderResult as Result } from '../result.js';
import { type SpzFile, SpzMagic, SpzShDims } from './schema.js';

async function parseInternal(ctx: RuntimeContext, input: ReadonlyUint8Array): Promise<Result<SpzFile>> {
  let data = input;
  const magicOf = (d: ReadonlyUint8Array) =>
    d.length >= 4 ? new DataView(d.buffer, d.byteOffset, d.byteLength).getUint32(0, true) : 0;
  // gzipped (versions 1 to 3), or inflated already (a download's)
  if (magicOf(data) !== SpzMagic && isGzip(data as Uint8Array)) {
    data = await ungzip(ctx, data as Uint8Array<ArrayBuffer>);
  }
  if (data.length < 16 || magicOf(data) !== SpzMagic) return Result.error('not an SPZ file (no NGSP header)');
  const view = new DataView(data.buffer, data.byteOffset, data.byteLength);
  const version = view.getUint32(4, true);
  // version 4: zstd streams, a header of its own
  if (version < 1 || version > 3) return Result.error(`SPZ version ${version} is not supported (versions 1 to 3)`);
  const count = view.getUint32(8, true);
  const shDegree = data[12];
  const dim = SpzShDims[shDegree];
  if (dim === undefined) return Result.error(`SPZ of spherical harmonics degree ${shDegree}`);

  let o = 16;
  const section = (bytes: number) => {
    const s = data.subarray(o, o + bytes);
    o += bytes;
    return s;
  };
  const positions = section(count * (version === 1 ? 6 : 9));
  const alphas = section(count);
  const colors = section(count * 3);
  const scales = section(count * 3);
  const rotations = section(count * (version >= 3 ? 4 : 3));
  const sh = section(count * dim * 3);
  if (o > data.length) return Result.error(`SPZ data ends before its ${count} splats`);
  return Result.success<SpzFile>({
    version,
    count,
    shDegree,
    fractionalBits: data[13],
    flags: data[14],
    positions,
    alphas,
    colors,
    scales,
    rotations,
    sh,
  });
}

/** Parses Niantic's `.spz` (versions 1 to 3, gzipped or not): its header, its arrays as views (no copy) */
export function parseSpz(data: ReadonlyUint8Array) {
  return Task.create<Result<SpzFile>>('Parse SPZ', (ctx) => parseInternal(ctx, data));
}
