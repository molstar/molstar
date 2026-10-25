/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Texture } from '@molstar/graphics/gl/texture';
import { createDataArray, type DataArray } from '@molstar/graphics/gl/renderable/util';
import { getGridDim, getTextureMeshSplatInstances } from '../color-smoothing.js';

// The instances the GPU smoothing of texture meshes draws (perf-audit 13/P1): the grid the shaders accumulate, emulated
// on the CPU (cells within 2 of a vertex, weighted `2 - distance`, summed in fp32 in draw order: instance-major, as
// instanced draws blend and the compute path's fixed-point sums give in any order), must be the same with the subset
// as with every instance.

type SplatInstances = Pick<
  ReturnType<typeof getTextureMeshSplatInstances>,
  'instanceCount' | 'transformBuffer' | 'instanceBuffer'
>;

const GroupCount = 7;
const VertexCount = 400;
const VertexBox = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(10, 10, 10));

/** Seeded pseudo-random numbers in [0, 1) */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 4294967296;
  };
}

function getVertices(withNegativeGroups: boolean) {
  const rnd = random(42);
  const positions = new Float32Array(VertexCount * 3);
  const groups = new Int32Array(VertexCount);
  for (let j = 0; j < VertexCount; ++j) {
    for (let k = 0; k < 3; ++k) positions[j * 3 + k] = rnd() * 10;
    // a cell whose voxel has no atom: group -1 (packed 0)
    groups[j] = withNegativeGroups && j % 50 === 0 ? -1 : Math.floor(rnd() * GroupCount);
  }
  return { positions, groups };
}

/** `count` instances along x, 12 apart (the cells of neighbours overlap), instance i at row 5i mod `count` */
function getInstances(count: number) {
  const transformBuffer = new Float32Array(count * 16);
  const instanceBuffer = new Float32Array(count);
  for (let i = 0; i < count; ++i) {
    Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(i * 12, 0, 0)), transformBuffer, i * 16);
    instanceBuffer[i] = (i * 5) % count;
  }
  return { instanceCount: count, transformBuffer, instanceBuffer };
}

function getInput(instances: SplatInstances, layer: DataArray<Uint8Array>, itemSize: 4 | 1) {
  const box = Box3D.setEmpty(Box3D());
  for (let i = 0; i < instances.instanceCount; ++i) {
    const t = Mat4.fromArray(Mat4(), instances.transformBuffer, i * 16);
    Box3D.addBox3D(box, Box3D.transform(Box3D(), VertexBox, t));
  }
  return {
    ...instances,
    vertexCount: VertexCount,
    groupCount: GroupCount,
    positionTexture: {} as Texture,
    groupTexture: {} as Texture,
    colorData: {} as Texture,
    colorType: 'groupInstance' as const,
    boundingSphere: Sphere3D.fromBox3D(Sphere3D(), box),
    invariantBoundingSphere: Sphere3D.fromBox3D(Sphere3D(), VertexBox),
    layer: { image: layer, itemSize },
    vertexBox: VertexBox,
  };
}

/** The grid as on the GPU paths: `Box3D.fromSphere3D` of the bounding sphere, expanded by `1 + resolution` */
function getGrid(input: ReturnType<typeof getInput>, resolution: number) {
  const box = Box3D.fromSphere3D(Box3D(), input.boundingSphere);
  const gridDim = getGridDim(box, resolution);
  const pad = 1 + resolution;
  const min = Vec3.sub(Vec3(), box.min, Vec3.create(pad, pad, pad));
  return { gridDim, min, scaleFactor: 1 / resolution };
}

/** The normalized grid of `instances` drawn (fp32 sums in draw order; no weight gives 0) */
function accumulate(
  instances: SplatInstances,
  vertices: ReturnType<typeof getVertices>,
  layer: DataArray<Uint8Array>,
  itemSize: number,
  grid: ReturnType<typeof getGrid>,
) {
  const { gridDim, min, scaleFactor } = grid;
  const [nx, ny, nz] = gridDim;
  const sums = new Float32Array(nx * ny * nz * itemSize);
  const weights = new Float32Array(nx * ny * nz);
  const t = Mat4();
  const v = Vec3();
  for (let i = 0; i < instances.instanceCount; ++i) {
    Mat4.fromArray(t, instances.transformBuffer, i * 16);
    const row = instances.instanceBuffer[i];
    for (let j = 0; j < VertexCount; ++j) {
      Vec3.fromArray(v, vertices.positions, j * 3);
      Vec3.transformMat4(v, v, t);
      Vec3.scale(v, Vec3.sub(v, v, min), scaleFactor);
      // an index before the layer (row 0, group -1) reads something
      const ci = (row * GroupCount + vertices.groups[j]) * itemSize;
      for (let z = Math.ceil(v[2] - 2); z <= Math.floor(v[2] + 2); ++z) {
        for (let y = Math.ceil(v[1] - 2); y <= Math.floor(v[1] + 2); ++y) {
          for (let x = Math.ceil(v[0] - 2); x <= Math.floor(v[0] + 2); ++x) {
            if (x < 0 || y < 0 || z < 0 || x >= nx || y >= ny || z >= nz) continue;
            const d = Math.hypot(x - v[0], y - v[1], z - v[2]);
            if (d > 2) continue;
            const f = 2 - d;
            const c = x + y * nx + z * nx * ny;
            for (let k = 0; k < itemSize; ++k) sums[c * itemSize + k] += (ci >= 0 ? layer.array[ci + k] : 255) * f;
            weights[c] += f;
          }
        }
      }
    }
  }
  const out = new Uint8Array(sums.length);
  for (let c = 0; c < weights.length; ++c) {
    if (weights[c] === 0) continue;
    for (let k = 0; k < itemSize; ++k) out[c * itemSize + k] = Math.round(sums[c * itemSize + k] / weights[c]);
  }
  return out;
}

/** Values that differ (not `toEqual`: its diff of large arrays runs out of memory) */
function countDiffs(a: Uint8Array, b: Uint8Array) {
  let n = a.length === b.length ? 0 : Infinity;
  for (let k = 0; k < a.length; ++k) if (a[k] !== b[k]) n += 1;
  return n;
}

/** A layer of `itemSize` that is non-zero in the rows `rows` (all groups), zero in the texels beyond */
function getLayer(instanceCount: number, itemSize: number, rows: number[]) {
  const layer = createDataArray(instanceCount * GroupCount, itemSize, Uint8Array);
  for (const r of rows) {
    for (let k = 0; k < GroupCount * itemSize; ++k) {
      layer.array[r * GroupCount * itemSize + k] = 50 + ((r * 31 + k) % 200);
    }
  }
  return layer;
}

describe('texture mesh smoothing instances', () => {
  it('draws the reached instances and their neighbours, in order, with the same grid', () => {
    const instances = getInstances(12);
    const vertices = getVertices(false);
    for (const itemSize of [4, 1] as const) {
      const layer = getLayer(12, itemSize, [6]);
      const input = getInput(instances, layer, itemSize);
      const grid = getGrid(input, 1);
      const subset = getTextureMeshSplatInstances(input, grid.gridDim, grid.min, grid.scaleFactor);

      // instances 0 (row 0, always drawn) and 1; 5, 6 (row 6) and 7; 10 and 11 (row 7, which vertices of group -1
      // would let read the last item of row 6)
      expect(Array.from(subset.instanceBuffer)).toEqual([0, 5, 1, 6, 11, 2, 7]);
      expect(subset.splatCount).toBe(7);
      expect(Array.from(subset.transformBuffer.subarray(16, 32))).toEqual(
        Array.from(instances.transformBuffer.subarray(16, 32)),
      );
      expect(Array.from(subset.transformBuffer.subarray(32, 48))).toEqual(
        Array.from(instances.transformBuffer.subarray(5 * 16, 6 * 16)),
      );
      const all = accumulate(instances, vertices, layer, itemSize, grid);
      expect(all.some((x) => x !== 0)).toBe(true);
      expect(countDiffs(accumulate(subset, vertices, layer, itemSize, grid), all)).toBe(0);
    }
  });

  it('counts the item before a row, which vertices of group -1 read', () => {
    const instances = getInstances(12);
    const vertices = getVertices(true);
    const layer = getLayer(12, 4, [6]);
    const input = getInput(instances, layer, 4);
    const grid = getGrid(input, 1);
    const subset = getTextureMeshSplatInstances(input, grid.gridDim, grid.min, grid.scaleFactor);

    // instance 11 (row 7), not a neighbour of instance 6 (row 6), reads the last item of row 6 through its vertices of
    // group -1
    expect(Array.from(subset.instanceBuffer)).toEqual([0, 5, 1, 6, 11, 2, 7]);
    expect(
      countDiffs(accumulate(subset, vertices, layer, 4, grid), accumulate(instances, vertices, layer, 4, grid)),
    ).toBe(0);
  });

  it('draws only row 0 and its neighbours for a zero layer', () => {
    const instances = getInstances(12);
    const vertices = getVertices(true);
    const layer = getLayer(12, 1, []);
    const input = getInput(instances, layer, 1);
    const grid = getGrid(input, 1);
    const subset = getTextureMeshSplatInstances(input, grid.gridDim, grid.min, grid.scaleFactor);

    // instance 0 (row 0, whose vertices of group -1 read before the layer) and 1
    expect(Array.from(subset.instanceBuffer)).toEqual([0, 5]);
    expect(
      countDiffs(accumulate(subset, vertices, layer, 1, grid), accumulate(instances, vertices, layer, 1, grid)),
    ).toBe(0);
  });

  it('draws every instance without a vertex box, a layer, or for a non-instance type', () => {
    const instances = getInstances(12);
    const input = getInput(instances, getLayer(12, 4, [6]), 4);
    const grid = getGrid(input, 1);
    const { gridDim, min, scaleFactor } = grid;
    for (const other of [
      { ...input, vertexBox: undefined },
      { ...input, layer: undefined },
      { ...input, colorType: 'group' as const },
    ]) {
      const subset = getTextureMeshSplatInstances(other, gridDim, min, scaleFactor);
      expect(subset.splatCount).toBe(12);
      expect(subset.transformBuffer).toBe(instances.transformBuffer);
      expect(subset.instanceBuffer).toBe(instances.instanceBuffer);
    }
  });
});
