/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { GridPositions } from '../grid-positions.js';

/** Positions over a box as large as the yeast cell's (about 47,000 Å), a fixed sequence */
function createPositions(n: number) {
  let s = 11;
  const rnd = () => (s = (s * 1664525 + 1013904223) >>> 0) / 4294967296;
  const out = new Float32Array(n * 3);
  for (let i = 0; i < n * 3; ++i) out[i] = (rnd() - 0.5) * 47000;
  out.set([-23500, -23500, -23500], 0);
  out.set([23500, 23500, 23500], 3);
  return out;
}

describe('GridPositions', () => {
  it('holds positions within half a step of the box (21 bits an axis), the box ends exact', () => {
    const positions = createPositions(1000);
    const grid = GridPositions.fromFloat(positions);
    expect(grid.data.length).toBe(2000);
    expect(GridPositions.count(grid)).toBe(1000);
    expect(grid.step).toBeCloseTo(47000 / GridPositions.Steps, 6);
    const a = Vec3();
    for (let i = 0; i < 1000; ++i) {
      GridPositions.get(grid, i, a);
      for (let k = 0; k < 3; ++k) {
        expect(Math.abs(a[k] - positions[i * 3 + k])).toBeLessThanOrEqual(grid.step / 2 + 1e-3 * grid.step);
      }
    }
    expect(Array.from(GridPositions.get(grid, 0, a))).toEqual([-23500, -23500, -23500]);
    // the yeast cell's extent: about 0.011 Å
    expect(grid.step / 2).toBeLessThan(0.012);
  });

  it('reads float32 positions as they are, and decodes a grid to float32', () => {
    const positions = createPositions(10);
    expect(GridPositions.is(positions)).toBe(false);
    expect(GridPositions.toFloat(positions)).toBe(positions);
    expect(Array.from(GridPositions.get(positions, 3, Vec3()))).toEqual(Array.from(positions.subarray(9, 12)));
    const grid = GridPositions.fromFloat(positions);
    expect(GridPositions.is(grid)).toBe(true);
    const decoded = GridPositions.toFloat(grid);
    for (let k = 0; k < decoded.length; ++k) expect(Math.abs(decoded[k] - positions[k])).toBeLessThan(grid.step);
  });

  it('writes a point and no positions', () => {
    const point = GridPositions.fromFloat(new Float32Array([1, 2, 3, 1, 2, 3]));
    expect(point.step).toBe(1);
    expect(Array.from(GridPositions.get(point, 1, Vec3()))).toEqual([1, 2, 3]);
    expect(GridPositions.count(GridPositions.fromFloat(new Float32Array(0)))).toBe(0);
  });
});
