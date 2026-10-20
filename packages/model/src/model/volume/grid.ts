/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type SpacegroupCell, Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { Tensor, Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import {
  type Histogram,
  type HistogramRobustStats,
  calculateHistogram,
  histogramRobustStats,
} from '@molstar/core/math/histogram';
import { lerp } from '@molstar/core/math/interpolate';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const v3lerp = Vec3.lerp;

/** The basic unit cell that contains the grid data. */
interface Grid {
  readonly transform: Grid.Transform;
  readonly cells: Tensor;
  readonly stats: Readonly<{
    min: number;
    max: number;
    mean: number;
    sigma: number;
  }>;
  readonly periodicity?: 'none' | 'xyz';
}

namespace Grid {
  export const One: Grid = {
    transform: { kind: 'matrix', matrix: Mat4.identity() },
    cells: Tensor.create(Tensor.Space([1, 1, 1], [0, 1, 2]), Tensor.Data1([0])),
    stats: { min: 0, max: 0, mean: 0, sigma: 0 },
  };

  export type Transform =
    | { kind: 'spacegroup'; cell: SpacegroupCell; fractionalBox: Box3D }
    | { kind: 'matrix'; matrix: Mat4 };

  const _scale = Mat4(),
    _translate = Mat4();
  export function getGridToCartesianTransform(grid: Grid) {
    if (grid.transform.kind === 'matrix') {
      return Mat4.copy(Mat4(), grid.transform.matrix);
    }

    if (grid.transform.kind === 'spacegroup') {
      const {
        cells: { space },
      } = grid;
      const scale = Mat4.fromScaling(
        _scale,
        Vec3.div(Vec3(), Box3D.size(Vec3(), grid.transform.fractionalBox), Vec3.ofArray(space.dimensions)),
      );
      const translate = Mat4.fromTranslation(_translate, grid.transform.fractionalBox.min);
      return Mat4.mul3(Mat4(), grid.transform.cell.fromFractional, translate, scale);
    }

    return Mat4.identity();
  }

  export function areEquivalent(gridA: Grid, gridB: Grid) {
    return gridA === gridB;
  }

  export function isEmpty(grid: Grid) {
    return grid.cells.data.length === 0;
  }

  export function getBoundingSphere(grid: Grid, boundingSphere?: Sphere3D) {
    if (!boundingSphere) boundingSphere = Sphere3D();

    const dimensions = grid.cells.space.dimensions as Vec3;
    const transform = Grid.getGridToCartesianTransform(grid);
    return Sphere3D.fromDimensionsAndTransform(boundingSphere, dimensions, transform);
  }

  const _isoBbox = Box3D();
  export function getIsosurfaceBoundingSphere(grid: Grid, isoValue: number, boundingSphere?: Sphere3D) {
    const neg = isoValue < 0;

    const c = [0, 0, 0];
    const getCoords = grid.cells.space.getCoords;
    const d = grid.cells.data;
    const [xn, yn, zn] = grid.cells.space.dimensions;

    let minx = xn - 1,
      miny = yn - 1,
      minz = zn - 1;
    let maxx = 0,
      maxy = 0,
      maxz = 0;
    for (let i = 0, il = d.length; i < il; ++i) {
      if ((neg && d[i] <= isoValue) || (!neg && d[i] >= isoValue)) {
        getCoords(i, c);
        if (c[0] < minx) minx = c[0];
        if (c[1] < miny) miny = c[1];
        if (c[2] < minz) minz = c[2];
        if (c[0] > maxx) maxx = c[0];
        if (c[1] > maxy) maxy = c[1];
        if (c[2] > maxz) maxz = c[2];
      }
    }

    Vec3.set(_isoBbox.min, minx - 1, miny - 1, minz - 1);
    Vec3.set(_isoBbox.max, maxx + 1, maxy + 1, maxz + 1);
    const transform = Grid.getGridToCartesianTransform(grid);
    Box3D.transform(_isoBbox, _isoBbox, transform);
    return Sphere3D.fromBox3D(boundingSphere || Sphere3D(), _isoBbox);
  }

  /**
   * Compute histogram with given bin count.
   * Cached on the Grid object.
   */
  export function getHistogram(grid: Grid, binCount: number) {
    let histograms = (grid as any)._histograms as { [binCount: number]: Histogram };
    if (!histograms) {
      histograms = (grid as any)._histograms = {};
    }
    if (!histograms[binCount]) {
      histograms[binCount] = calculateHistogram(grid.cells.data, binCount, {
        min: grid.stats.min,
        max: grid.stats.max,
      });
    }
    return histograms[binCount];
  }

  /** Default bin count for histogram-derived robust statistics. */
  export const RobustStatsBinCount = 1024;

  /**
   * Compute robust summary statistics (percentiles, mean, sigma) from the
   * cached histogram. Set `ignoreZero` to skip the histogram bin that
   * contains exact zero values (typical for masked/padded density maps).
   * Cached on the Grid object.
   */
  export function getRobustStats(grid: Grid, opts?: { ignoreZero?: boolean; binCount?: number }): HistogramRobustStats {
    const binCount = opts?.binCount ?? RobustStatsBinCount;
    const ignoreZero = !!opts?.ignoreZero;
    const key = `${binCount}|${ignoreZero ? 1 : 0}`;
    let cache = (grid as any)._robustStats as { [key: string]: HistogramRobustStats };
    if (!cache) cache = (grid as any)._robustStats = {};
    if (!cache[key]) {
      cache[key] = histogramRobustStats(getHistogram(grid, binCount), { ignoreZero });
    }
    return cache[key];
  }

  export function makeGetTrilinearlyInterpolated(grid: Grid, transform: 'none' | 'relative') {
    const cartnToGrid = Grid.getGridToCartesianTransform(grid);
    Mat4.invert(cartnToGrid, cartnToGrid);
    const gridCoords = Vec3();

    const { stats } = grid;
    const { dimensions, get } = grid.cells.space;
    const data = grid.cells.data;

    const [mi, mj, mk] = dimensions;
    const getValue = (i: number, j: number, k: number) => get(data, i, j, k);
    const isPeriodic = grid.periodicity === 'xyz';

    return function getTrilinearlyInterpolated(position: Vec3): number {
      v3transformMat4(gridCoords, position, cartnToGrid);

      const value = trilinearlyInterpolate(gridCoords, mi, mj, mk, isPeriodic, getValue);
      if (Number.isNaN(value)) return value;

      if (transform === 'relative') {
        return (value - stats.mean) / stats.sigma;
      } else {
        return value;
      }
    };
  }

  /**
   * Core trilinear interpolation function.
   * @param gridCoords - Position in grid coordinates (fractional indices)
   * @param mi, mj, mk - Grid dimensions
   * @param getValue - Function to get value at integer grid coordinates
   * @returns Interpolated value or NaN if out of bounds
   */
  export function trilinearlyInterpolate(
    gridCoords: Vec3,
    mi: number,
    mj: number,
    mk: number,
    isPeriodic: boolean,
    getValue: (i: number, j: number, k: number) => number,
  ): number {
    const i = Math.trunc(gridCoords[0]);
    const j = Math.trunc(gridCoords[1]);
    const k = Math.trunc(gridCoords[2]);

    if (!isPeriodic && (i < 0 || i >= mi || j < 0 || j >= mj || k < 0 || k >= mk)) {
      return Number.NaN;
    }

    const u = gridCoords[0] - i;
    const v = gridCoords[1] - j;
    const w = gridCoords[2] - k;

    const ii = Math.min(i + 1, mi - 1);
    const jj = Math.min(j + 1, mj - 1);
    const kk = Math.min(k + 1, mk - 1);

    let a = getValue(i, j, k);
    let b = getValue(ii, j, k);
    let c = getValue(i, jj, k);
    let d = getValue(ii, jj, k);
    const x = lerp(lerp(a, b, u), lerp(c, d, u), v);

    a = getValue(i, j, kk);
    b = getValue(ii, j, kk);
    c = getValue(i, jj, kk);
    d = getValue(ii, jj, kk);
    const y = lerp(lerp(a, b, u), lerp(c, d, u), v);

    return lerp(x, y, w);
  }

  const _a = Vec3(),
    _b = Vec3(),
    _c = Vec3(),
    _d = Vec3();
  const _ab = Vec3(),
    _cd = Vec3(),
    _x = Vec3(),
    _y = Vec3();

  /**
   * Core trilinear interpolation function for Vec3 values.
   * More efficient for interleaved data like gradients since getValue is called once per grid point.
   * @param gridCoords - Position in grid coordinates (fractional indices)
   * @param mi, mj, mk - Grid dimensions
   * @param getValue - Function to get Vec3 value at integer grid coordinates
   * @param out - Output Vec3
   * @returns true if interpolation succeeded, false if out of bounds
   */
  export function trilinearlyInterpolateVec3(
    gridCoords: Vec3,
    mi: number,
    mj: number,
    mk: number,
    getValue: (i: number, j: number, k: number, out: Vec3) => void,
    out: Vec3,
  ): boolean {
    const i = Math.trunc(gridCoords[0]);
    const j = Math.trunc(gridCoords[1]);
    const k = Math.trunc(gridCoords[2]);

    if (i < 0 || i >= mi || j < 0 || j >= mj || k < 0 || k >= mk) {
      return false;
    }

    const u = gridCoords[0] - i;
    const v = gridCoords[1] - j;
    const w = gridCoords[2] - k;

    const ii = Math.min(i + 1, mi - 1);
    const jj = Math.min(j + 1, mj - 1);
    const kk = Math.min(k + 1, mk - 1);

    // Interpolate in the k plane
    getValue(i, j, k, _a);
    getValue(ii, j, k, _b);
    getValue(i, jj, k, _c);
    getValue(ii, jj, k, _d);
    v3lerp(_ab, _a, _b, u);
    v3lerp(_cd, _c, _d, u);
    v3lerp(_x, _ab, _cd, v);

    // Interpolate in the k+1 plane
    getValue(i, j, kk, _a);
    getValue(ii, j, kk, _b);
    getValue(i, jj, kk, _c);
    getValue(ii, jj, kk, _d);
    v3lerp(_ab, _a, _b, u);
    v3lerp(_cd, _c, _d, u);
    v3lerp(_y, _ab, _cd, v);

    // Final interpolation between planes
    v3lerp(out, _x, _y, w);
    return true;
  }

  type Gradients = {
    values: Float32Array;
    magnitude: {
      min: number;
      max: number;
      mean: number;
      sigma: number;
    };
  };

  /**
   * The gradient at a voxel by central differences (one-sided at the faces, divisor 1 on an axis of size 1). Reads the
   * data directly where the space allows it (the same values as `get`), else through `get`.
   */
  function makeGetGradient(grid: Grid): (i: number, j: number, k: number, out: Vec3) => void {
    const { space, data } = grid.cells;
    const [mi, mj, mk] = space.dimensions;
    const strides = Tensor.strides(space);
    if (strides && strides.length === 3) {
      const [sx, sy, sz] = strides;
      return (i, j, k, out) => {
        const im = Math.max(0, i - 1);
        const ip = Math.min(mi - 1, i + 1);
        const jm = Math.max(0, j - 1);
        const jp = Math.min(mj - 1, j + 1);
        const km = Math.max(0, k - 1);
        const kp = Math.min(mk - 1, k + 1);
        const o = i * sx + j * sy + k * sz;
        out[0] = (data[o + (ip - i) * sx] - data[o + (im - i) * sx]) / (ip - im || 1);
        out[1] = (data[o + (jp - j) * sy] - data[o + (jm - j) * sy]) / (jp - jm || 1);
        out[2] = (data[o + (kp - k) * sz] - data[o + (km - k) * sz]) / (kp - km || 1);
      };
    }
    const { get } = space;
    return (i, j, k, out) => {
      const im = Math.max(0, i - 1);
      const ip = Math.min(mi - 1, i + 1);
      const jm = Math.max(0, j - 1);
      const jp = Math.min(mj - 1, j + 1);
      const km = Math.max(0, k - 1);
      const kp = Math.min(mk - 1, k + 1);
      out[0] = (get(data, ip, j, k) - get(data, im, j, k)) / (ip - im || 1);
      out[1] = (get(data, i, jp, k) - get(data, i, jm, k)) / (jp - jm || 1);
      out[2] = (get(data, i, j, kp) - get(data, i, j, km)) / (kp - km || 1);
    };
  }

  /**
   * Gradients at every voxel (`makeGetGradient`), laid out as the data: xyz interleaved at `dataOffset * 3`, with the
   * magnitudes' statistics. Not cached: 12 bytes a voxel, the caller keeps it as long as it needs it.
   */
  export function getGradients(grid: Grid): Gradients {
    const { dimensions, dataOffset } = grid.cells.space;
    const [mi, mj, mk] = dimensions;
    const getGradient = makeGetGradient(grid);
    const n = mi * mj * mk;
    const values = new Float32Array(n * 3);
    const g = Vec3();

    let min = Infinity;
    let max = -Infinity;
    let sum = 0;
    let sumSq = 0;

    for (let k = 0; k < mk; ++k) {
      for (let j = 0; j < mj; ++j) {
        for (let i = 0; i < mi; ++i) {
          getGradient(i, j, k, g);
          const idx = dataOffset(i, j, k) * 3;
          values[idx] = g[0];
          values[idx + 1] = g[1];
          values[idx + 2] = g[2];

          const mag = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
          if (mag < min) min = mag;
          if (mag > max) max = mag;
          sum += mag;
          sumSq += mag * mag;
        }
      }
    }

    if (min === Infinity) min = 0;
    if (max === -Infinity) max = 1;
    return {
      values,
      magnitude: { min: Math.sqrt(min), max: Math.sqrt(max), mean: sum / n, sigma: Math.sqrt(sumSq / n) },
    };
  }

  /**
   * Create a function that returns the trilinearly interpolated gradient at a grid position. The gradients of the eight
   * surrounding voxels are computed on demand, nothing is kept: a gradient field cached for the whole grid was 12 bytes
   * a voxel for the volume's lifetime (perf-audit 04/P8). They are rounded to float32 as that field stored them.
   */
  export function makeGetInterpolatedGradient(grid: Grid) {
    const [mi, mj, mk] = grid.cells.space.dimensions;
    const getGradient = makeGetGradient(grid);
    const fround = Math.fround;
    const getGradientVec3 = (i: number, j: number, k: number, out: Vec3) => {
      getGradient(i, j, k, out);
      out[0] = fround(out[0]);
      out[1] = fround(out[1]);
      out[2] = fround(out[2]);
    };

    return function getInterpolatedGradient(gridCoords: Vec3, out: Vec3): boolean {
      return trilinearlyInterpolateVec3(gridCoords, mi, mj, mk, getGradientVec3, out);
    };
  }
}

export { Grid };
