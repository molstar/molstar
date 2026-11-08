/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';
import { SegmentedArray } from './segmented-array.js';

/**
 * Affine transforms on the GPU: an attribute of item size 12 holds per item the first three rows of a column-major
 * `Mat4` (16 floats on the CPU) as three vec4s (`aTransform0` to `aTransform2`); the shaders rebuild the matrix with the
 * last row 0, 0, 0, 1 (GLSL `affineFromRows`, WGSL `affineTransform`)
 */
export const AffineRowsItemSize = 12;

/**
 * Packs the column-major `Mat4`s of `src` (16 floats each) as their first three rows (12 floats each) into `out`, from
 * `outOffset`, an item each `outStride` floats (12: packed)
 */
export function packAffineRows(
  src: ArrayLike<number>,
  count: number,
  out: Float32Array,
  outOffset = 0,
  outStride = AffineRowsItemSize,
) {
  // unrolled: millions of instances (mesoscale)
  for (let i = 0, s = 0, o = outOffset; i < count; ++i, s += 16, o += outStride) {
    out[o] = src[s];
    out[o + 1] = src[s + 4];
    out[o + 2] = src[s + 8];
    out[o + 3] = src[s + 12];
    out[o + 4] = src[s + 1];
    out[o + 5] = src[s + 5];
    out[o + 6] = src[s + 9];
    out[o + 7] = src[s + 13];
    out[o + 8] = src[s + 2];
    out[o + 9] = src[s + 6];
    out[o + 10] = src[s + 10];
    out[o + 11] = src[s + 14];
  }
  return out;
}

/**
 * `packAffineRows` of the first `count` transforms of `src`, one array or a merged object's members' (`SegmentedArray`,
 * each part packed at its offset: not composed into one array first, millions of instances)
 */
export function packAffineRowsOf(
  src: Float32Array | SegmentedArray<Float32Array>,
  count: number,
  out: Float32Array,
  outOffset = 0,
  outStride = AffineRowsItemSize,
) {
  if (!SegmentedArray.is(src)) return packAffineRows(src, count, out, outOffset, outStride);
  for (const p of src.parts) {
    const first = p.offset / 16;
    const n = Math.min(p.count / 16, count - first);
    if (n > 0) packAffineRows(p.array, n, out, outOffset + first * outStride, outStride);
  }
  return out;
}

/**
 * A cell with the transforms of `cell` as rows (`packAffineRows`), repacked when its version changes; values other than
 * a `Float32Array` or a `SegmentedArray` of them (instanced attributes written on the GPU) pass through; observers are
 * those of `cell`
 */
export function affineRowsCell<T>(cell: ValueCell<T>): ValueCell<T> {
  const ref = { id: cell.ref.id, version: -1, value: undefined as unknown };
  return {
    get observers() {
      return cell.observers;
    },
    set observers(value) {
      cell.observers = value;
    },
    get ref() {
      const source = cell.ref;
      if (source.version !== ref.version || source.id !== ref.id) {
        const v = source.value as unknown;
        if (v instanceof Float32Array || (SegmentedArray.is(v) && v.ctor === Float32Array)) {
          const count = Math.floor(v.length / 16);
          const out =
            ref.value instanceof Float32Array && ref.value.length === count * 12
              ? ref.value
              : new Float32Array(count * 12);
          ref.value = packAffineRowsOf(v as Float32Array | SegmentedArray<Float32Array>, count, out);
        } else {
          ref.value = v;
        }
        ref.id = source.id;
        ref.version = source.version;
      }
      return ref;
    },
  } as unknown as ValueCell<T>;
}
