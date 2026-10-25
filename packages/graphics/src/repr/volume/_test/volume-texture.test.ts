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

/** the texture as built cell by cell through `dataOffset` with `Vec3` (the implementation before strides) */
function reference(volume: Volume, type: 'byte' | 'float' | 'halfFloat') {
  const {
    cells: { space, data },
    stats: { max, min },
  } = volume.grid;
  const [width, height, depth] = space.dimensions as Vec3;
  const o = space.dataOffset;
  const array =
    type === 'byte'
      ? new Uint8Array(width * height * depth * 4)
      : type === 'halfFloat'
        ? new Uint16Array(width * height * depth * 4)
        : new Float32Array(width * height * depth * 4);
  const diff = max - min;
  const n0 = Vec3();
  const n1 = Vec3();
  let i = 0;
  for (let z = 0; z < depth; ++z) {
    for (let y = 0; y < height; ++y) {
      for (let x = 0; x < width; ++x) {
        const offset = o(x, y, z);
        Vec3.set(
          n0,
          data[o(Math.max(0, x - 1), y, z)],
          data[o(x, Math.max(0, y - 1), z)],
          data[o(x, y, Math.max(0, z - 1))],
        );
        Vec3.set(
          n1,
          data[o(Math.min(width - 1, x + 1), y, z)],
          data[o(x, Math.min(height - 1, y + 1), z)],
          data[o(x, y, Math.min(depth - 1, z + 1))],
        );
        Vec3.normalize(n0, Vec3.sub(n0, n0, n1));
        Vec3.addScalar(n0, Vec3.scale(n0, n0, 0.5), 0.5);
        if (type === 'byte') {
          Vec3.toArray(Vec3.scale(n0, n0, 255), array, i);
          array[i + 3] = Math.round(((data[offset] - min) / diff) * 255);
        } else if (type === 'halfFloat') {
          array[i] = toHalfFloat(n0[0]);
          array[i + 1] = toHalfFloat(n0[1]);
          array[i + 2] = toHalfFloat(n0[2]);
          array[i + 3] = toHalfFloat((data[offset] - min) / diff);
        } else {
          Vec3.toArray(n0, array, i);
          array[i + 3] = (data[offset] - min) / diff;
        }
        i += 4;
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
