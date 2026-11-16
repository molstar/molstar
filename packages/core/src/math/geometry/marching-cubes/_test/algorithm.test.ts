/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Tensor } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { computeMarchingCubes, getMarchingCubesInputParams, type MarchingCubesBuilder } from '../algorithm.js';
import { estimateIsosurfaceVertexCount } from '../estimate.js';

/** A density falling off from the grid center, `n` minus the distance: a sphere of radius `r` at iso value `n - r` */
function sphereField(n: number) {
  const values = new Float32Array(n * n * n);
  const c = (n - 1) / 2;
  for (let z = 0; z < n; ++z) {
    for (let y = 0; y < n; ++y) {
      for (let x = 0; x < n; ++x) values[x + y * n + z * n * n] = n - Math.hypot(x - c, y - c, z - c);
    }
  }
  // x fastest
  return { values, scalarField: Tensor.create(Tensor.Space([n, n, n], [2, 1, 0], Float32Array), Tensor.Data1(values)) };
}

/** Plain arrays: what a builder of core gets, without a geometry */
function arraysBuilder(): MarchingCubesBuilder<{
  vertices: number[];
  normals: number[];
  groups: number[];
  triangles: number[];
}> {
  const vertices: number[] = [];
  const normals: number[] = [];
  const groups: number[] = [];
  const triangles: number[] = [];
  return {
    addVertex: (x, y, z) => {
      vertices.push(x, y, z);
      return vertices.length / 3 - 1;
    },
    addNormal: (x, y, z) => {
      normals.push(x, y, z);
    },
    addGroup: (g) => {
      groups.push(g);
    },
    addTriangle: (vertList, a, b, c) => {
      const i = vertList[a];
      const j = vertList[b];
      const k = vertList[c];
      if (i >= 0 && j >= 0 && k >= 0) triangles.push(i, j, k);
    },
    get: () => ({ vertices, normals, groups, triangles }),
  };
}

describe('marching cubes', () => {
  it('makes a closed sphere through any builder: vertices on the iso level, every edge shared by two triangles', async () => {
    const n = 12;
    const radius = 4;
    const { values, scalarField } = sphereField(n);
    const params = getMarchingCubesInputParams({ isoLevel: n - radius, scalarField });
    const { vertices, normals, triangles } = await computeMarchingCubes(
      RuntimeContext.Synchronous,
      arraysBuilder(),
      params,
    );

    const vertexCount = vertices.length / 3;
    expect(vertexCount).toBeGreaterThan(0);
    expect(normals.length).toBe(vertices.length);
    // the vertices interpolated on the cells' edges: about the radius from the center
    const c = (n - 1) / 2;
    for (let i = 0; i < vertexCount; ++i) {
      const d = Math.hypot(vertices[i * 3] - c, vertices[i * 3 + 1] - c, vertices[i * 3 + 2] - c);
      expect(Math.abs(d - radius)).toBeLessThan(0.1);
    }
    // closed: each undirected edge in exactly two triangles
    const edges = new Map<string, number>();
    for (let t = 0; t < triangles.length; t += 3) {
      for (const [a, b] of [
        [triangles[t], triangles[t + 1]],
        [triangles[t + 1], triangles[t + 2]],
        [triangles[t + 2], triangles[t]],
      ]) {
        const key = a < b ? `${a}-${b}` : `${b}-${a}`;
        edges.set(key, (edges.get(key) ?? 0) + 1);
      }
    }
    expect([...edges.values()].every((count) => count === 2)).toBe(true);
    // as the estimate of the compute extraction counts them (three per triangle, exact at this size)
    expect(triangles.length).toBe(estimateIsosurfaceVertexCount(values, [n, n, n], [1, n, n * n], n - radius));
  });
});
