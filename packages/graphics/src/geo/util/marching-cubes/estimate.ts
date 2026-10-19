/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { TriTable } from './tables.js';

/** Triangles of each cube configuration (`TriTable`) */
const TriangleCounts = Uint8Array.from(TriTable, (t) => t.length / 3);

/** Cells visited at most (about 6 ms); larger grids are sampled by z-slices */
const MaxCells = 2 ** 20;
/** Headroom of a sampled estimate over the count it scales */
const SampledHeadroom = 1.2;

/**
 * Vertex count of the marching cubes isosurface of `data` at `isoLevel`, unindexed (three per triangle, as the compute
 * extraction writes, see `computeChunkedIsosurface`): the triangles of the cells (`TriTable`, corners at or below
 * `isoLevel` as the compute kernels' `cubeCase`). Exact for grids up to `MaxCells` cells; larger ones visit every n-th slice of cells
 * along z, scaled, with headroom. `dim` is the grid's point count per axis, `strides` the data offsets of a step along
 * each axis, `dataDim` the data's point count per axis where the grid wraps around it (periodic volumes, the grid one
 * point larger, as the compute kernels' `value`).
 */
export function estimateIsosurfaceVertexCount(
  data: ArrayLike<number>,
  dim: ArrayLike<number>,
  strides: ArrayLike<number>,
  isoLevel: number,
  dataDim: ArrayLike<number> = dim,
): number {
  const [nx, ny, nz] = [dim[0], dim[1], dim[2]];
  if (nx < 2 || ny < 2 || nz < 2) return 0;
  // the data offset of each grid position along an axis
  const offsets = (n: number, d: number, stride: number) => {
    const o = new Array<number>(n);
    for (let i = 0; i < n; ++i) o[i] = (i < d ? i : i - d) * stride;
    return o;
  };
  const [ox, oy, oz] = [
    offsets(nx, dataDim[0], strides[0]),
    offsets(ny, dataDim[1], strides[1]),
    offsets(nz, dataDim[2], strides[2]),
  ];
  const cellSlices = nz - 1;
  const sliceCells = (nx - 1) * (ny - 1);
  const step = Math.max(1, Math.ceil((sliceCells * cellSlices) / MaxCells));

  // per point of a slice: whether it is at or below `isoLevel`
  const lower = new Uint8Array(nx * ny);
  const upper = new Uint8Array(nx * ny);
  const below = (out: Uint8Array, k: number) => {
    for (let j = 0; j < ny; ++j) {
      const row = j * nx;
      const offset = oy[j] + oz[k];
      for (let i = 0; i < nx; ++i) out[row + i] = data[offset + ox[i]] <= isoLevel ? 1 : 0;
    }
  };

  let triangles = 0;
  let sampled = 0;
  for (let k = 0; k < cellSlices; k += step) {
    below(lower, k);
    below(upper, k + 1);
    for (let j = 0; j < ny - 1; ++j) {
      const r0 = j * nx;
      const r1 = r0 + nx;
      for (let i = 0; i < nx - 1; ++i) {
        // the corner order of `TriTable` (`cubeCase`'s `Corners`)
        const index =
          lower[r0 + i] |
          (lower[r0 + i + 1] << 1) |
          (lower[r1 + i + 1] << 2) |
          (lower[r1 + i] << 3) |
          (upper[r0 + i] << 4) |
          (upper[r0 + i + 1] << 5) |
          (upper[r1 + i + 1] << 6) |
          (upper[r1 + i] << 7);
        triangles += TriangleCounts[index];
      }
    }
    sampled += 1;
  }
  const vertices = 3 * triangles;
  return step === 1 ? vertices : Math.ceil(((vertices * cellSlices) / sampled) * SampledHeadroom);
}
