/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4 } from '../3d/mat4.js';
import { MinimizeRmsd, RmsdTransformState } from '../3d/minimize-rmsd.js';
import { Vec3 } from '../3d/vec3.js';

function positions(xyz: number[][]): MinimizeRmsd.Positions {
  return { x: xyz.map((p) => p[0]), y: xyz.map((p) => p[1]), z: xyz.map((p) => p[2]) };
}

const a1 = positions([
  [0, 0, 0],
  [1, 0, 0],
  [0, 2, 0],
  [0, 0, 3],
]);
const b1 = positions([
  [5, 0, 0],
  [6, 0, 0],
  [5, 2, 0],
  [5, 0, 3],
]);
const a2 = positions([
  [10, 10, 10],
  [11, 10, 10],
  [10, 12, 10],
  [10, 10, 13],
]);
const b2 = positions([
  [-3, 1, 0],
  [-3, 2, 0],
  [-5, 1, 0],
  [-3, 1, 3],
]);

describe('MinimizeRmsd', () => {
  it('a reused state neither writes into nor keeps a caller centre (perf-audit-2 02/B4)', () => {
    const centerA = Vec3.create(0.25, 0.5, 0.75);
    const centerB = Vec3.create(5.25, 0.5, 0.75);
    const first = { a: a1, b: b1, centerA, centerB };
    const state = new RmsdTransformState(first, { bTransform: Mat4(), rmsd: 0, nAlignedElements: 0 });
    MinimizeRmsd.compute(first, undefined, state);

    const reused = MinimizeRmsd.compute({ a: a2, b: b2 }, undefined, state);
    expect(Vec3.equals(centerA, Vec3.create(0.25, 0.5, 0.75))).toBe(true);
    expect(Vec3.equals(centerB, Vec3.create(5.25, 0.5, 0.75))).toBe(true);

    const fresh = MinimizeRmsd.compute({ a: a2, b: b2 });
    expect(reused.rmsd).toBeCloseTo(fresh.rmsd, 10);
    expect(Mat4.areEqual(reused.bTransform, fresh.bTransform, 1e-9)).toBe(true);
  });

  it('gives the input its computed centres', () => {
    const input: MinimizeRmsd.Input = { a: a1, b: b1 };
    const r = MinimizeRmsd.compute(input);
    expect(r.rmsd).toBeCloseTo(0, 6);
    expect(Vec3.equals(input.centerA!, Vec3.create(0.25, 0.5, 0.75))).toBe(true);
  });
});
