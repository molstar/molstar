/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { getSceneGrid, SceneCellClasses } from '../compute/scene-cells.js';

const dims = (grid: ReturnType<typeof getSceneGrid>, k: number) => [1, 2, 3].map((a) => grid.grids[k * 4 + a]);

describe('getSceneGrid', () => {
  it('spans the members bounds, a class per size, cells for the instances per cell', () => {
    const bounds = [
      Sphere3D.create(Vec3.create(0, 0, 0), 50),
      Sphere3D.create(Vec3.create(100, 0, 0), 50),
      // no instances: neither its bounds nor its class count
      Sphere3D.create(Vec3.create(1000, 0, 0), 50),
    ];
    const grid = getSceneGrid(bounds, [0.5, 20, 400], [1000, 10, 0], 1010, 10);
    expect(Array.from(grid.min.slice(0, 3))).toEqual([-50, -50, -50]);
    // the box 200 x 100 x 100 Å, 101 cells of 10 instances: a side of about 27 Å
    const side = 1 / grid.grids[0];
    expect(side).toBeCloseTo(Math.cbrt((200 * 100 * 100) / 101), 6);
    expect(dims(grid, 0)).toEqual([8, 4, 4]);
    // a diameter of 40 Å: the next class (cells of twice the side)
    expect(grid.classes[1]).toBe(1);
    expect(1 / grid.grids[4]).toBeCloseTo(side * 2, 6);
    expect(dims(grid, 1)).toEqual([4, 2, 2]);
    // past the largest class's side the largest; a class without instances no cells
    expect(grid.classes[2]).toBe(SceneCellClasses - 1);
    expect(dims(grid, 2)).toEqual([0, 0, 0]);
    expect(dims(grid, 3)).toEqual([0, 0, 0]);
    expect(Array.from(grid.firsts)).toEqual([0, 128, 144, 144]);
    expect(grid.cellCount).toBe(144);
  });

  it('grows the cells of a class past its most cells; no cells without instances', () => {
    const grid = getSceneGrid([Sphere3D.create(Vec3(), 1e4)], [0], [1e9], 1e9, 1);
    const [nx, ny, nz] = dims(grid, 0);
    expect(nx * ny * nz).toBeLessThanOrEqual(1 << 22);
    expect(nx * ny * nz).toBeGreaterThan((1 << 22) * 0.9);

    const empty = getSceneGrid([], [], [], 0, 32);
    expect(empty.cellCount).toBe(0);
    expect(Array.from(empty.min)).toEqual([0, 0, 0, 0]);
  });
});
