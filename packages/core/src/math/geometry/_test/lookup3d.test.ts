/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { GridLookup3D } from '../../geometry.js';
import { sortArray } from '@molstar/core/data/util/sort';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { getBoundary } from '../boundary.js';

const xs = [0, 0, 1];
const ys = [0, 1, 0];
const zs = [0, 0, 0];
const rs = [0, 0.5, 1 / 3];

describe('GridLookup3d', () => {
  it('basic', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    let r = grid.find(0, 0, 0, 0);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(0);

    r = grid.nearest(0, 0, 0, 1);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(0);

    r = grid.find(0, 0, 0, 1);
    expect(r.count).toBe(3);
    expect(sortArray(r.indices)).toEqual([0, 1, 2]);

    r = grid.nearest(0, 0, 0, 3);
    expect(r.count).toBe(3);
    expect(sortArray(r.indices)).toEqual([0, 1, 2]);
  });

  it('radius', () => {
    const position = { x: xs, y: ys, z: zs, radius: [0, 0.5, 1 / 3], indices: OrderedSet.ofBounds(0, 3) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    let r = grid.find(0, 0, 0, 0);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(0);

    r = grid.nearest(0, 0, 0, 1);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(0);

    r = grid.find(0, 0, 0, 0.5);
    expect(r.count).toBe(2);
    expect(sortArray(r.indices)).toEqual([0, 1]);

    r = grid.nearest(0, 0, 0, 3);
    expect(r.count).toBe(3);
    expect(sortArray(r.indices)).toEqual([0, 1, 2]);
  });

  it('indexed', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofSingleton(1), radius: rs };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    let r = grid.find(0, 0, 0, 0);
    expect(r.count).toBe(0);

    r = grid.nearest(0, 0, 0, 1);
    expect(r.count).toBe(1);

    r = grid.find(0, 0, 0, 0.5);
    expect(r.count).toBe(1);
    expect(sortArray(r.indices)).toEqual([0]);

    r = grid.nearest(0, 0, 0, 3);
    expect(r.count).toBe(1);
    expect(sortArray(r.indices)).toEqual([0]);
  });

  it('sorted-array indices', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofSortedArray([0, 2]), radius: rs };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    let r = grid.find(0, 0, 0, 0);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(0);

    // element 1 (y = 1) is excluded, so within radius 1 only elements 0 and 2 are found
    r = grid.find(0, 0, 0, 1);
    expect(r.count).toBe(2);
    expect(sortArray(r.indices)).toEqual([0, 1]);

    r = grid.nearest(0, 0, 0, 2);
    expect(r.count).toBe(2);
    expect(sortArray(r.indices)).toEqual([0, 1]);
  });

  it('check', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    expect(grid.check(0, 0, 0, 0.5)).toBe(true);
    expect(grid.check(5, 5, 5, 0.5)).toBe(false);
  });

  it('approxNearest', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    let r = grid.approxNearest(0, 0, 0, 0.5);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(0);

    r = grid.approxNearest(5, 5, 5, 0.5);
    expect(r.count).toBe(0);
  });

  it('approxNearest gives the nearest of a cell, once (perf-audit-2 02/B1)', () => {
    // one cell: the farther element first in the bucket, then the nearer one
    const position = { x: [2.9, 0.3, 4], y: [0, 0, 0], z: [0, 0, 0], indices: OrderedSet.ofBounds(0, 3) };
    const grid = GridLookup3D(position, getBoundary(position), Vec3.create(10, 10, 10));
    const r = grid.approxNearest(0, 0, 0, 5);
    expect(r.count).toBe(1);
    expect(r.indices[0]).toBe(1);
    expect(r.squaredDistances[0]).toBeCloseTo(0.09, 10);
  });

  it('empty', () => {
    const position = { x: [] as number[], y: [] as number[], z: [] as number[], indices: OrderedSet.ofBounds(0, 0) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    expect(grid.find(0, 0, 0, 1).count).toBe(0);
    expect(grid.nearest(0, 0, 0, 1).count).toBe(0);
    expect(grid.check(0, 0, 0, 1)).toBe(false);
  });

  it('nearest from outside the boundary', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    const r = grid.nearest(50, 50, 50, 1);
    expect(r.count).toBe(1);
    // element 1 (0,1,0) and element 2 (1,0,0) are equidistant candidates; (0,1,0) wins on distance ties by heap order
    expect(r.squaredDistances[0]).toBeCloseTo(
      Math.min((50 - 0) ** 2 + (50 - 1) ** 2 + 50 ** 2, (50 - 1) ** 2 + 50 ** 2 + 50 ** 2),
    );
  });

  it('nearest with stopIf', () => {
    const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
    const boundary = getBoundary(position);
    const grid = GridLookup3D(position, boundary);

    const r = grid.nearest(0, 0, 0, 3, () => true);
    expect(r.count).toBeGreaterThanOrEqual(1);
    expect(r.count).toBeLessThanOrEqual(3);
  });

  it('many points with fixed cell size match brute force', () => {
    // deterministic LCG so the test is reproducible
    let seed = 42;
    const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;

    const count = 1000;
    const px = new Float32Array(count);
    const py = new Float32Array(count);
    const pz = new Float32Array(count);
    for (let i = 0; i < count; ++i) {
      px[i] = rand() * 100;
      py[i] = rand() * 100;
      pz[i] = rand() * 100;
    }

    const position = { x: px, y: py, z: pz, indices: OrderedSet.ofBounds(0, count) };
    const boundary = getBoundary(position);
    // sparse fixed-size cells exercise the occupied-cell compaction path
    const grid = GridLookup3D(position, boundary, Vec3.create(3, 3, 3));

    const queries: [number, number, number, number][] = [
      [50, 50, 50, 10],
      [0, 0, 0, 25],
      [100, 100, 100, 15],
      [50, 0, 100, 5],
    ];
    for (const [qx, qy, qz, qr] of queries) {
      const expected: number[] = [];
      for (let i = 0; i < count; ++i) {
        const dx = px[i] - qx,
          dy = py[i] - qy,
          dz = pz[i] - qz;
        if (dx * dx + dy * dy + dz * dz <= qr * qr) expected.push(i);
      }
      const r = grid.find(qx, qy, qz, qr);
      expect(sortArray(r.indices.slice(0, r.count))).toEqual(expected);
      expect(grid.check(qx, qy, qz, qr)).toBe(expected.length > 0);
    }
  });
});

// Pins the exact hit order and squared distances: intra-unit bonds, contacts and the molecular-surface torus loop
// consume hits in the order the grid returns them, so a layout change must not reorder or round them.
describe('GridLookup3d result order', () => {
  function points(count: number, seed: number) {
    const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;
    const x = new Float64Array(count);
    const y = new Float64Array(count);
    const z = new Float64Array(count);
    const r = new Float64Array(count);
    for (let i = 0; i < count; ++i) {
      x[i] = rand() * 30;
      y[i] = rand() * 20;
      z[i] = rand() * 25;
      r[i] = 1 + rand();
    }
    return { x, y, z, r };
  }

  // FNV-1a over (index, float64 bits of the squared distance) of every hit, in order
  const f64 = new Float64Array(1);
  const u32 = new Uint32Array(f64.buffer);
  function fold(h: number, v: number) {
    return Math.imul(h ^ v, 16777619) >>> 0;
  }
  function foldResult(
    h: number,
    r: { count: number; indices: ArrayLike<number>; squaredDistances: ArrayLike<number> },
  ) {
    h = fold(h, r.count);
    for (let i = 0; i < r.count; i++) {
      f64[0] = r.squaredDistances[i];
      h = fold(fold(fold(h, r.indices[i]), u32[0]), u32[1]);
    }
    return h;
  }

  const count = 2000;
  const { x, y, z, r } = points(count, 7);
  // every third element, so lookup indices differ from the element indices
  const sorted = OrderedSet.ofSortedArray(Int32Array.from({ length: Math.ceil(count / 3) }, (_, i) => 3 * i));
  const cases = {
    'f64 interval': { x, y, z, indices: OrderedSet.ofBounds(0, count) },
    'f32 interval': {
      x: Float32Array.from(x),
      y: Float32Array.from(y),
      z: Float32Array.from(z),
      indices: OrderedSet.ofBounds(0, count),
    },
    'array sorted': { x: Array.from(x), y: Array.from(y), z: Array.from(z), indices: sorted },
    'f64 interval radii': { x, y, z, radius: r, indices: OrderedSet.ofBounds(0, count) },
    'f32 sorted radii': {
      x: Float32Array.from(x),
      y: Float32Array.from(y),
      z: Float32Array.from(z),
      radius: Float32Array.from(r),
      indices: sorted,
    },
  };

  function fingerprint(position: (typeof cases)[keyof typeof cases], cell?: Vec3) {
    const grid = GridLookup3D(position, getBoundary(position), cell);
    let h = 0x811c9dc5;
    let hits = 0;
    for (let i = 0; i < 300; i++) {
      // query points on a lattice that covers the box and its surroundings
      const qx = -3 + (i % 7) * 6.1;
      const qy = -3 + (Math.floor(i / 7) % 6) * 4.7;
      const qz = -3 + Math.floor(i / 42) * 4.3;
      for (const qr of [0.8, 2, 4, 6.4]) {
        const f = grid.find(qx, qy, qz, qr);
        hits += f.count;
        h = foldResult(h, f);
        h = foldResult(h, grid.approxNearest(qx, qy, qz, qr));
        h = fold(h, +grid.check(qx, qy, qz, qr));
      }
      h = foldResult(h, grid.nearest(qx, qy, qz, 1));
      h = foldResult(h, grid.nearest(qx, qy, qz, 6));
    }
    return `${hits} ${h.toString(16)}`;
  }

  it('find hits keep their order and distances', () => {
    const position = {
      x: x.subarray(0, 40),
      y: y.subarray(0, 40),
      z: z.subarray(0, 40),
      indices: OrderedSet.ofBounds(0, 40),
    };
    const grid = GridLookup3D(position, getBoundary(position));
    const f = grid.find(15, 10, 12, 8);
    expect([
      Array.from(f.indices.slice(0, f.count)),
      Array.from(f.squaredDistances.slice(0, f.count)),
    ]).toMatchInlineSnapshot(`
      [
        [
          29,
          39,
          14,
          27,
          12,
          0,
          13,
          15,
        ],
        [
          63.09301965765512,
          49.87555500551978,
          48.0678555329588,
          54.76454431533104,
          28.444327327750038,
          27.02092338930176,
          52.993755782136034,
          8.150830517972967,
        ],
      ]
    `);
  });

  it('find, approxNearest, nearest and check fingerprints', () => {
    const out: Record<string, string> = {};
    for (const [name, position] of Object.entries(cases)) {
      out[name] = fingerprint(position);
      out[`${name}, 3 Å cells`] = fingerprint(position, Vec3.create(3, 3, 3));
    }
    expect(out).toMatchInlineSnapshot(`
      {
        "array sorted": "7261 18c73e1e",
        "array sorted, 3 Å cells": "7261 e8f28869",
        "f32 interval": "21810 4c64db64",
        "f32 interval, 3 Å cells": "21810 ba8b3fd",
        "f32 sorted radii": "15459 f64e7abc",
        "f32 sorted radii, 3 Å cells": "15459 2105f01e",
        "f64 interval": "21810 551a81f3",
        "f64 interval radii": "46167 954d18fa",
        "f64 interval radii, 3 Å cells": "46167 137f36a",
        "f64 interval, 3 Å cells": "21810 6c89a522",
      }
    `);
  });

  it('a contiguous SortedArray gives the same results as the interval', () => {
    const contiguous = OrderedSet.ofSortedArray(Int32Array.from({ length: count }, (_, i) => i));
    expect(fingerprint({ ...cases['f32 interval'], indices: contiguous })).toBe(fingerprint(cases['f32 interval']));
    expect(fingerprint({ ...cases['f64 interval radii'], indices: contiguous })).toBe(
      fingerprint(cases['f64 interval radii']),
    );
    // offset range: elements 100..1099 as a SortedArray vs as an interval
    const sub = OrderedSet.ofSortedArray(Int32Array.from({ length: 1000 }, (_, i) => 100 + i));
    expect(fingerprint({ ...cases['f32 interval'], indices: sub })).toBe(
      fingerprint({ ...cases['f32 interval'], indices: OrderedSet.ofBounds(100, 1100) }),
    );
  });
});
