/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { toHalfFloat } from '@molstar/core/util/number-conversion';
import { CustomProperties } from '@molstar/model/model/custom-property';
import type { Volume } from '@molstar/model/model/volume';
import { createVolumeTexture3d } from '../util.js';

/** a volume of `dim` cells with pseudo-random values, flat (0) where x is below 2 */
function createVolume(dim: [number, number, number], axisOrder: [number, number, number]): Volume {
  const space = Tensor.Space(dim, axisOrder, Float32Array);
  const data = new Float32Array(dim[0] * dim[1] * dim[2]);
  let seed = 7;
  for (let x = 0; x < dim[0]; ++x) {
    for (let y = 0; y < dim[1]; ++y) {
      for (let z = 0; z < dim[2]; ++z) {
        seed = (seed * 16807) % 2147483647;
        data[space.dataOffset(x, y, z)] = x < 2 ? 0 : (seed / 2147483647) * 4 - 1;
      }
    }
  }
  let min = Infinity;
  let max = -Infinity;
  for (const v of data) {
    min = Math.min(min, v);
    max = Math.max(max, v);
  }
  return {
    label: 'test',
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, Tensor.Data1(data)),
      stats: { min, max, mean: 0, sigma: 0 },
    },
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'test', name: 'test', data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

/** the texture as built cell by cell through `dataOffset`: the values normalized by the grid's min and max, x fastest */
function reference(volume: Volume, type: 'byte' | 'float' | 'halfFloat') {
  const {
    cells: { space, data },
    stats: { max, min },
  } = volume.grid;
  const [width, height, depth] = space.dimensions as Vec3;
  const o = space.dataOffset;
  const count = width * height * depth;
  const array =
    type === 'byte' ? new Uint8Array(count) : type === 'halfFloat' ? new Uint16Array(count) : new Float32Array(count);
  const diff = max - min;
  let i = 0;
  for (let z = 0; z < depth; ++z) {
    for (let y = 0; y < height; ++y) {
      for (let x = 0; x < width; ++x, ++i) {
        const value = (data[o(x, y, z)] - min) / diff;
        array[i] = type === 'byte' ? Math.round(value * 255) : type === 'halfFloat' ? toHalfFloat(value) : value;
      }
    }
  }
  return array;
}

describe('createVolumeTexture3d', () => {
  for (const type of ['byte', 'float', 'halfFloat'] as const) {
    it(`matches the cell by cell build (${type}, every axis order)`, () => {
      for (const axisOrder of [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
      ] as [number, number, number][]) {
        const volume = createVolume([7, 5, 6], axisOrder);
        expect(Array.from(createVolumeTexture3d(volume, type).array)).toEqual(Array.from(reference(volume, type)));
      }
    });
  }

  it('reads spaces that override their accessors through them (no strides: periodic, wrapped volumes)', () => {
    const volume = createVolume([7, 5, 6], [2, 1, 0]);
    const { space } = volume.grid.cells;
    const nx = space.dimensions[0];
    // x mirrored: not linear in x
    const mirrored = { ...space, dataOffset: (x: number, y: number, z: number) => space.dataOffset(nx - 1 - x, y, z) };
    const v = { ...volume, grid: { ...volume.grid, cells: { ...volume.grid.cells, space: mirrored } } } as Volume;
    expect(Tensor.strides(mirrored)).toBeUndefined();
    expect(Array.from(createVolumeTexture3d(v, 'float').array)).toEqual(Array.from(reference(v, 'float')));
  });

  it('handles grids one cell thick', () => {
    for (const dim of [
      [1, 4, 3],
      [4, 1, 3],
      [4, 3, 1],
      [1, 1, 1],
    ] as [number, number, number][]) {
      const volume = createVolume(dim, [2, 1, 0]);
      expect(Array.from(createVolumeTexture3d(volume, 'float').array)).toEqual(Array.from(reference(volume, 'float')));
    }
  });
});
