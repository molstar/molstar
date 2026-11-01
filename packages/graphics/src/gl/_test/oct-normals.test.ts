/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { octNormalsCell, packOctNormals } from '../oct-normals.js';
import { SegmentedArray } from '../segmented-array.js';

/** The shaders' `octDecode`, normalized as its callers do */
function decode(a: number, b: number) {
  const u = a / 32767;
  const v = b / 32767;
  const z = 1 - Math.abs(u) - Math.abs(v);
  const t = Math.max(-z, 0);
  const x = u + (u >= 0 ? -t : t);
  const y = v + (v >= 0 ? -t : t);
  const l = Math.hypot(x, y, z);
  return [x / l, y / l, z / l];
}

describe('oct normals', () => {
  it('encode unit vectors within 1e-4 radians, both hemispheres and the axes', () => {
    const normals: number[] = [1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1];
    let seed = 1;
    const random = () => ((seed = (seed * 16807) % 2147483647) / 2147483647) * 2 - 1;
    for (let i = 0; i < 1000; ++i) {
      const [x, y, z] = [random(), random(), random()];
      const l = Math.hypot(x, y, z);
      normals.push(x / l, y / l, z / l);
    }
    const count = normals.length / 3;
    const encoded = packOctNormals(normals, count, new Int16Array(count * 2));
    let worst = 0;
    for (let i = 0; i < count; ++i) {
      const [x, y, z] = decode(encoded[i * 2], encoded[i * 2 + 1]);
      const dot = x * normals[i * 3] + y * normals[i * 3 + 1] + z * normals[i * 3 + 2];
      worst = Math.max(worst, Math.acos(Math.min(1, dot)));
    }
    expect(worst).toBeLessThan(1e-4);
  });

  it('keep a segmented array encoded part by part, its identity while the layout stays', () => {
    const a = new Float32Array([0, 0, 1, 0, 0, 1]);
    const b = new Float32Array([1, 0, 0]);
    const parts = [
      { array: a, offset: 0, count: 6, bias: 0, version: 0 },
      { array: b, offset: 6, count: 3, bias: 0, version: 0 },
    ];
    const cell = ValueCell.create(SegmentedArray.create(Float32Array, parts));
    const encoded = octNormalsCell(cell);
    const first = encoded.ref.value as unknown as SegmentedArray<Int16Array>;
    expect(first.length).toBe(6);
    expect(first.parts.map((p) => [p.offset, p.count])).toEqual([
      [0, 4],
      [4, 2],
    ]);
    const versions = first.parts.map((p) => p.version);

    // the second member's normals change: the same array, only its part encoded again
    parts[1].array = new Float32Array([0, 1, 0]);
    parts[1].version += 1;
    ValueCell.update(cell, cell.ref.value);
    const second = encoded.ref.value as unknown as SegmentedArray<Int16Array>;
    expect(second).toBe(first);
    expect(second.parts[0].version).toBe(versions[0]);
    expect(second.parts[1].version).toBe(versions[1] + 1);
    expect(Array.from(second.parts[1].array)).toEqual([0, 32767]);
  });
});
