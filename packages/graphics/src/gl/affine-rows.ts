/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';

/**
 * Affine transforms on the GPU: an attribute of item size 12 holds per item the first three rows of a column-major
 * `Mat4` (16 floats on the CPU) as three vec4s (`aTransform0` to `aTransform2`); the shaders rebuild the matrix with the
 * last row 0, 0, 0, 1 (GLSL `affineFromRows`, WGSL `affineTransform`)
 */
export const AffineRowsItemSize = 12;

/** Packs the column-major `Mat4`s of `src` (16 floats each) as their first three rows (12 floats each) into `out` */
export function packAffineRows(src: ArrayLike<number>, count: number, out: Float32Array, outOffset = 0) {
  for (let i = 0; i < count; ++i) {
    const s = i * 16;
    const o = outOffset + i * 12;
    for (let r = 0; r < 3; ++r) {
      out[o + r * 4] = src[s + r];
      out[o + r * 4 + 1] = src[s + 4 + r];
      out[o + r * 4 + 2] = src[s + 8 + r];
      out[o + r * 4 + 3] = src[s + 12 + r];
    }
  }
  return out;
}

/**
 * A cell with the transforms of `cell` as rows (`packAffineRows`), repacked when its version changes; values other than
 * a `Float32Array` (instanced attributes written on the GPU) pass through; observers are those of `cell`
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
        if (v instanceof Float32Array) {
          const count = Math.floor(v.length / 16);
          const out =
            ref.value instanceof Float32Array && ref.value.length === count * 12
              ? ref.value
              : new Float32Array(count * 12);
          ref.value = packAffineRows(v, count, out);
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
