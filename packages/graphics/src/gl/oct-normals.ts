/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ValueCell } from '@molstar/core/util';
import { SegmentedArray } from './segmented-array.js';

/**
 * Unit vectors on the GPU in two signed normalized 16-bit integers (octahedral encoding, `oct16` attributes: 4 instead
 * of 12 bytes a vertex): on the CPU three floats (`Mesh.normalBuffer`); the shaders decode them (GLSL and WGSL
 * `octDecode`, then normalize), within about 1e-4 radians
 */
export const OctNormalItemSize = 3;

const signNotZero = (v: number) => (v >= 0 ? 1 : -1);

/** Encodes the `count` unit vectors of `src` (three floats each) into `out` (two `Int16` each, at `outOffset`) */
export function packOctNormals(src: ArrayLike<number>, count: number, out: Int16Array, outOffset = 0) {
  for (let i = 0; i < count; ++i) {
    const x = src[i * 3];
    const y = src[i * 3 + 1];
    const z = src[i * 3 + 2];
    const l1 = Math.abs(x) + Math.abs(y) + Math.abs(z);
    let u = l1 > 0 ? x / l1 : 0;
    let v = l1 > 0 ? y / l1 : 0;
    if (z < 0) {
      const fu = (1 - Math.abs(v)) * signNotZero(u);
      v = (1 - Math.abs(u)) * signNotZero(v);
      u = fu;
    }
    out[outOffset + i * 2] = Math.round(Math.max(-1, Math.min(1, u)) * 32767);
    out[outOffset + i * 2 + 1] = Math.round(Math.max(-1, Math.min(1, v)) * 32767);
  }
  return out;
}

/** `packOctNormals` of all of `src` into a new array */
function encode(src: Float32Array) {
  const count = Math.floor(src.length / 3);
  return packOctNormals(src, count, new Int16Array(count * 2));
}

/**
 * A cell with the normals of `cell` encoded (`packOctNormals`), again when its version changes. A segmented array (the
 * normals of merged meshes, their members' arrays as parts) maps to one of encoded parts that keeps its identity while
 * the parts' layout does, a part encoded again where its version changed, so that uploads write only those
 * (`SegmentedArray.writeChanged`); observers are those of `cell`
 */
export function octNormalsCell<T>(cell: ValueCell<T>): ValueCell<T> {
  const ref = { id: cell.ref.id, version: -1, value: undefined as unknown };
  let source: SegmentedArray | undefined;
  let encoded: SegmentedArray<Int16Array> | undefined;
  /** per part the source array and version it was encoded from */
  let encodedFrom: { array: unknown; version: number }[] = [];
  const encodeSegmented = (s: SegmentedArray) => {
    if (s !== source || !encoded) {
      source = s;
      encodedFrom = s.parts.map(() => ({ array: undefined, version: -1 }));
      const parts = s.parts.map((p) => ({
        array: new Int16Array(0),
        offset: (p.offset / 3) * 2,
        count: (p.count / 3) * 2,
        bias: 0,
        version: 0,
      }));
      encoded = SegmentedArray.create(Int16Array, parts, (s.length / 3) * 2);
    }
    for (let i = 0, il = s.parts.length; i < il; ++i) {
      const p = s.parts[i];
      const from = encodedFrom[i];
      if (from.array === p.array && from.version === p.version) continue;
      const e = encoded.parts[i];
      e.array = packOctNormals(p.array, p.count / 3, new Int16Array((p.count / 3) * 2));
      e.version += 1;
      from.array = p.array;
      from.version = p.version;
    }
    return encoded;
  };
  return {
    get observers() {
      return cell.observers;
    },
    set observers(value) {
      cell.observers = value;
    },
    get ref() {
      const box = cell.ref;
      if (box.version !== ref.version || box.id !== ref.id) {
        const v = box.value as unknown;
        if (v instanceof Float32Array) ref.value = encode(v);
        else if (SegmentedArray.is(v)) ref.value = encodeSegmented(v);
        else ref.value = v;
        ref.id = box.id;
        ref.version = box.version;
      }
      return ref;
    },
  } as unknown as ValueCell<T>;
}
