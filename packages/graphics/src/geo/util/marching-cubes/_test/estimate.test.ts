/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Tensor } from '@molstar/core/math/linear-algebra';
import { computeMarchingCubesMesh } from '../algorithm.js';
import { estimateIsosurfaceVertexCount } from '../estimate.js';
import { createWrappedTensor } from '../../../../repr/volume/util.js';

/** A density falling off from the grid center, `n` minus the distance: a sphere of radius `r` at iso value `n - r` */
function sphereValues(n: number) {
  const values = new Float32Array(n * n * n);
  const c = (n - 1) / 2;
  for (let z = 0; z < n; ++z) {
    for (let y = 0; y < n; ++y) {
      for (let x = 0; x < n; ++x) values[x + y * n + z * n * n] = n - Math.hypot(x - c, y - c, z - c);
    }
  }
  return values;
}

/** A periodic density, `n` minus the distance to `c` across the faces: a sphere of radius `r` at iso value `n - r`, cut by
 * the faces near `c` and continued on the opposite ones */
function periodicSphereValues(n: number, c: number) {
  const values = new Float32Array(n * n * n);
  const d = (v: number) => Math.min(Math.abs(v - c), n - Math.abs(v - c));
  for (let z = 0; z < n; ++z) {
    for (let y = 0; y < n; ++y) {
      for (let x = 0; x < n; ++x) values[x + y * n + z * n * n] = n - Math.hypot(d(x), d(y), d(z));
    }
  }
  return values;
}

/** Unindexed vertices of the CPU marching cubes (three per triangle, as the compute extraction writes) */
async function cpuVertexCount(values: Float32Array, n: number, iso: number) {
  // x fastest, as the strides below
  const scalarField = Tensor.create(Tensor.Space([n, n, n], [2, 1, 0], Float32Array), Tensor.Data1(values));
  const mesh = await computeMarchingCubesMesh({ isoLevel: iso, scalarField }).run();
  return mesh.triangleCount * 3;
}

describe('estimateIsosurfaceVertexCount', () => {
  it('counts the vertices of a grid it visits whole exactly', async () => {
    const n = 32;
    const values = sphereValues(n);
    // between grid values, where at or below and below agree
    const iso = n - 10.25;
    const estimate = estimateIsosurfaceVertexCount(values, [n, n, n], [1, n, n * n], iso);
    expect(estimate).toBe(await cpuVertexCount(values, n, iso));
  });

  it('follows strided data', async () => {
    const n = 24;
    const values = sphereValues(n);
    // the same grid with x slowest (z fastest): transposed data, transposed strides
    const transposed = new Float32Array(values.length);
    for (let z = 0; z < n; ++z) {
      for (let y = 0; y < n; ++y) {
        for (let x = 0; x < n; ++x) transposed[z + y * n + x * n * n] = values[x + y * n + z * n * n];
      }
    }
    const iso = n - 7.25;
    expect(estimateIsosurfaceVertexCount(transposed, [n, n, n], [n * n, n, 1], iso)).toBe(
      estimateIsosurfaceVertexCount(values, [n, n, n], [1, n, n * n], iso),
    );
  });

  it('samples large grids by slices, with headroom over the count', async () => {
    // 119^3 cells, more than visited whole
    const n = 120;
    const values = sphereValues(n);
    const iso = n - 40.25;
    const estimate = estimateIsosurfaceVertexCount(values, [n, n, n], [1, n, n * n], iso);
    const count = await cpuVertexCount(values, n, iso);
    expect(estimate).toBeGreaterThanOrEqual(count);
    expect(estimate).toBeLessThan(count * 1.3);
  });

  it('wraps a grid one point larger around the data', async () => {
    const n = 20;
    const values = periodicSphereValues(n, 1.5);
    const iso = n - 6.25;
    // the wrapped tensor of periodic volumes, as the CPU isosurface extracts it
    const tensor = Tensor.create(Tensor.Space([n, n, n], [2, 1, 0], Float32Array), Tensor.Data1(values));
    const mesh = await computeMarchingCubesMesh({ isoLevel: iso, scalarField: createWrappedTensor(tensor) }).run();
    const m = n + 1;
    const estimate = estimateIsosurfaceVertexCount(values, [m, m, m], [1, n, n * n], iso, [n, n, n]);
    expect(estimate).toBe(mesh.triangleCount * 3);
  });

  it('is zero without cells or crossings', () => {
    expect(estimateIsosurfaceVertexCount(new Float32Array(4), [4, 1, 1], [1, 4, 4], 0.5)).toBe(0);
    expect(estimateIsosurfaceVertexCount(new Float32Array(27), [3, 3, 3], [1, 3, 9], 0.5)).toBe(0);
  });
});
