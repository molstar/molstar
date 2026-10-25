/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor } from '@molstar/core/math/linear-algebra';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { Volume } from '@molstar/model/model/volume';
import { createVolumeBrickRange } from '../util.js';

/** a volume of `dim` cells, all 0 but the cells given as [x, y, z, value]; its min 0 and max 1 */
function createVolume(
  dim: [number, number, number],
  cells: [number, number, number, number][],
  axisOrder: [number, number, number] = [2, 1, 0],
): Volume {
  const space = Tensor.Space(dim, axisOrder, Float32Array);
  const data = new Float32Array(dim[0] * dim[1] * dim[2]);
  for (const [x, y, z, v] of cells) data[space.dataOffset(x, y, z)] = v;
  return {
    label: 'test',
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, Tensor.Data1(data)),
      stats: { min: 0, max: 1, mean: 0, sigma: 0 },
    },
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

describe('createVolumeBrickRange', () => {
  it('has a brick per 8 cells along each axis, the last one partial', () => {
    const { array, dim } = createVolumeBrickRange(createVolume([17, 8, 9], []));
    expect(Array.from(dim)).toEqual([3, 1, 2]);
    expect(array.length).toBe(3 * 1 * 2 * 2);
  });

  it('ranges over a brick cells, rounded outwards', () => {
    const { array } = createVolumeBrickRange(createVolume([8, 8, 8], [[3, 4, 5, 0.5]]));
    // 0.5 * 255 = 127.5
    expect(Array.from(array)).toEqual([0, 128]);
  });

  it('includes the next cell along each axis, which trilinear samples in the brick read', () => {
    // x = 8 is the first cell of the second brick along x and the next cell of the first
    const { array } = createVolumeBrickRange(createVolume([16, 8, 8], [[8, 0, 0, 1]]));
    expect(Array.from(array)).toEqual([0, 255, 0, 255]);
    // x = 9 is only in the second brick
    const second = createVolumeBrickRange(createVolume([16, 8, 8], [[9, 0, 0, 1]]));
    expect(Array.from(second.array)).toEqual([0, 1, 0, 255]);
  });

  it('is the same for any axis order of the data', () => {
    const cells: [number, number, number, number][] = [
      [8, 3, 1, 1],
      [2, 9, 16, 0.25],
      [17, 16, 8, 0.75],
    ];
    const expected = Array.from(createVolumeBrickRange(createVolume([18, 17, 19], cells)).array);
    for (const order of [
      [0, 1, 2],
      [1, 2, 0],
      [2, 0, 1],
    ] as [number, number, number][]) {
      expect(Array.from(createVolumeBrickRange(createVolume([18, 17, 19], cells, order)).array)).toEqual(expected);
    }
  });

  it('reads spaces that override their accessors through them (no strides: periodic, wrapped volumes)', () => {
    const cells: [number, number, number, number][] = [
      [3, 2, 1, 1],
      [12, 9, 4, 0.5],
    ];
    const dim: [number, number, number] = [18, 10, 9];
    const volume = createVolume(dim, cells);
    const { space } = volume.grid.cells;
    // x mirrored: not linear in x
    const mirrored = {
      ...space,
      dataOffset: (x: number, y: number, z: number) => space.dataOffset(dim[0] - 1 - x, y, z),
    };
    const v = { ...volume, grid: { ...volume.grid, cells: { ...volume.grid.cells, space: mirrored } } } as Volume;
    const expected = createVolume(
      dim,
      cells.map(([x, y, z, value]) => [dim[0] - 1 - x, y, z, value] as [number, number, number, number]),
    );
    expect(Array.from(createVolumeBrickRange(v).array)).toEqual(Array.from(createVolumeBrickRange(expected).array));
  });

  it('shares a corner cell with the seven bricks before it', () => {
    const { array } = createVolumeBrickRange(createVolume([16, 16, 16], [[8, 8, 8, 1]]));
    for (let i = 0; i < 8; ++i) expect(array[2 * i + 1]).toBe(255);
  });
});
