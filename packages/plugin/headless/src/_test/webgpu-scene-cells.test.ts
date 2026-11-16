/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util/value-cell';
import type { ComputeKernel, StorageBuffer } from '@molstar/graphics/gl/compute';
import { getDispatchSize } from '@molstar/graphics/gl/compute/prefix-sum';
import {
  createSceneCellBuilder,
  getSceneGrid,
  SceneCellKernels,
  type SceneGrid,
} from '@molstar/graphics/gl/compute/scene-cells';
import { CulledData } from '@molstar/graphics/gl/culled-data';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

/** A small deterministic generator (an LCG), values in [0, 1) */
function random(seed: number) {
  let s = seed >>> 0;
  return () => {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    return s / 2 ** 32;
  };
}

/**
 * Two segments: matrices (scaled unevenly, a small invariant sphere) and poses (rotated, scaled, a large one: another
 * size class), spread over a box shifted by `offset`
 */
function createScene(offset: Vec3) {
  const rnd = random(7);
  const n0 = 500;
  const n1 = 300;
  const matrices = new Float32Array(n0 * 16);
  const m = Mat4();
  for (let i = 0; i < n0; ++i) {
    Mat4.fromScaling(m, Vec3.create(0.5 + rnd(), 1, 0.5 + rnd() * 2));
    Mat4.setTranslation(m, Vec3.create(offset[0] + rnd() * 80, offset[1] + rnd() * 60, offset[2] + rnd() * 40));
    Mat4.toArray(m, matrices, i * 16);
  }
  const positions = new Float32Array(n1 * 3);
  const rotations = new Float32Array(n1 * 4);
  const scales = new Float32Array(n1);
  const q = Quat();
  for (let e = 0; e < n1; ++e) {
    positions.set([offset[0] + rnd() * 80, offset[1] + rnd() * 60, offset[2] + rnd() * 40], e * 3);
    Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(rnd() - 0.5, rnd(), 1)), rnd() * 6);
    Quat.toArray(q, rotations, e * 4);
    scales[e] = 0.8 + rnd() * 0.4;
  }
  const invariants = [Sphere3D.create(Vec3.create(0.5, 0, -0.25), 1), Sphere3D.create(Vec3.create(0, 1, 0), 30)];
  const transforms = [
    InstanceTransforms.ofMatrices(matrices, n0),
    InstanceTransforms.ofPoses({ positions, rotations, scales, count: n1, first: 0, center: invariants[1].center }),
  ];
  const segments = { count: 2, instanceBases: [0, n0], instanceCounts: [n0, n1] };
  return { transforms, invariants, segments, n: n0 + n1 };
}

type SceneInput = ReturnType<typeof createScene>;

/** each instance's sphere (f64, from its affine rows) and segment */
function getInstanceSpheres({ transforms, invariants, segments, n }: SceneInput) {
  const spheres = new Float64Array(n * 4);
  const segmentOf = new Uint32Array(n);
  const rows = new Float32Array(12);
  for (let s = 0; s < segments.count; ++s) {
    const { center: c, radius } = invariants[s];
    for (let i = 0; i < segments.instanceCounts[s]; ++i) {
      InstanceTransforms.writeAffineRows(transforms[s], i, 1, rows);
      const o = (segments.instanceBases[s] + i) * 4;
      for (let a = 0; a < 3; ++a) {
        spheres[o + a] = rows[a * 4] * c[0] + rows[a * 4 + 1] * c[1] + rows[a * 4 + 2] * c[2] + rows[a * 4 + 3];
      }
      let scale = 0;
      for (let col = 0; col < 3; ++col) {
        scale = Math.max(scale, Math.hypot(rows[col], rows[4 + col], rows[8 + col]));
      }
      spheres[o + 3] = radius * scale;
      segmentOf[segments.instanceBases[s] + i] = s;
    }
  }
  return { spheres, segmentOf };
}

/** the box of the spheres as a sphere and each segment's largest radius, as a merged object's bounds give them */
function getBounds(input: SceneInput, spheres: Float64Array) {
  const bounds: Sphere3D[] = [];
  const radii = new Float32Array(input.segments.count);
  for (let s = 0; s < input.segments.count; ++s) {
    const lo = [Infinity, Infinity, Infinity];
    const hi = [-Infinity, -Infinity, -Infinity];
    const base = input.segments.instanceBases[s];
    for (let i = 0; i < input.segments.instanceCounts[s]; ++i) {
      const o = (base + i) * 4;
      for (let a = 0; a < 3; ++a) {
        lo[a] = Math.min(lo[a], spheres[o + a] - spheres[o + 3]);
        hi[a] = Math.max(hi[a], spheres[o + a] + spheres[o + 3]);
      }
      radii[s] = Math.max(radii[s], spheres[o + 3]);
    }
    const center = Vec3.create((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2);
    bounds.push(Sphere3D.create(center, Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2));
  }
  return { bounds, radii };
}

/** the cell of a center in its class's grid, and whether it lies within `eps` of a cell's side */
function getCell(grid: SceneGrid, k: number, x: number, y: number, z: number, eps: number) {
  const inv = grid.grids[k * 4];
  const dims = [grid.grids[k * 4 + 1], grid.grids[k * 4 + 2], grid.grids[k * 4 + 3]];
  const g = [x, y, z].map((v, a) => (v - grid.min[a]) * inv);
  const near = g.some((v) => Math.abs(v - Math.round(v)) < eps);
  const c = g.map((v, a) => Math.min(dims[a] - 1, Math.max(0, Math.floor(v))));
  return { cell: grid.firsts[k] + c[0] + dims[0] * (c[1] + dims[1] * c[2]), near };
}

describe.skipIf(!WebGPU.module)('scene cells built on the GPU', () => {
  it('bins each instance into its class grid cell, each cell the sphere of its instances', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const compute = graphics.compute!;
      await compute.prepareKernels(SceneCellKernels);
      const builder = createSceneCellBuilder(compute);
      const cell = (bytes: number) => ValueCell.create<StorageBuffer>(compute.createBuffer(Math.max(16, bytes)));
      const dispatchItems = (kernel: ComputeKernel, n: number) => {
        // two rows of workgroups: the kernels linearize the id
        const { x, y } = getDispatchSize(Math.ceil(n / 64), 7);
        kernel.dispatch(x, y);
      };

      // built twice into the same buffers, the second after the instances moved: the boxes cleared in between
      for (const offset of [Vec3.create(-10, 5, 2), Vec3.create(120, -40, 30)]) {
        const input = createScene(offset);
        const { n, segments } = input;
        const { spheres, segmentOf } = getInstanceSpheres(input);
        const { bounds, radii } = getBounds(input, spheres);
        const grid = getSceneGrid(bounds, radii, segments.instanceCounts, n, 2);
        // the large invariant sphere a class of its own
        expect(grid.classes[1]).toBeGreaterThan(grid.classes[0]);

        const parts = input.transforms.map((transforms, s) => ({ transforms, first: segments.instanceBases[s] }));
        const layout = CulledData.layout(parts, new Uint32Array(n), segments, undefined, n, 128);
        const data = cell(layout.vec4s * 16);
        layout.write((d, o) => data.ref.value.write(d, o));
        const segmentData = cell(segments.count * 48);
        const segmentValues = new Float32Array(segments.count * 12);
        for (let s = 0; s < segments.count; ++s) {
          const { center, radius } = input.invariants[s];
          segmentValues.set([center[0], center[1], center[2], radius], s * 12);
        }
        segmentData.ref.value.write(segmentValues);
        const instanceCells = cell(n * 8);
        const cellSpheres = cell(grid.cellCount * 16);

        builder.set(grid, { instanceCount: n, segmentData, data, instanceCells, cellSpheres });
        builder.clear();
        compute.beginComputePass({ label: 'scene cells' });
        builder.dispatch(dispatchItems);
        compute.endComputePass();
        const cells = new Uint32Array(await instanceCells.ref.value.readAsync(0, n * 8));
        const gpuSpheres = new Float32Array(await cellSpheres.ref.value.readAsync(0, grid.cellCount * 16));

        // each instance in the cell of its center, of its segment's class (one off only on a cell's side, rounding)
        const members = new Map<number, number[]>();
        let nearSide = 0;
        for (let i = 0; i < n; ++i) {
          const s = segmentOf[i];
          expect(cells[i * 2 + 1]).toBe(s);
          const o = i * 4;
          const { cell: expected, near } = getCell(
            grid,
            grid.classes[s],
            spheres[o],
            spheres[o + 1],
            spheres[o + 2],
            1e-4,
          );
          if (near) nearSide += 1;
          else expect(cells[i * 2]).toBe(expected);
          const d = cells[i * 2];
          if (!members.has(d)) members.set(d, []);
          members.get(d)!.push(i);
        }
        expect(nearSide).toBeLessThan(n / 50);
        expect(members.size).toBeGreaterThan(20);

        // each cell's sphere that of the box of its instances' centers, by their largest radius; the others empty
        for (let d = 0; d < grid.cellCount; ++d) {
          const r = gpuSpheres[d * 4 + 3];
          const instances = members.get(d);
          if (!instances) {
            expect(r).toBeLessThan(0);
            continue;
          }
          const lo = [Infinity, Infinity, Infinity];
          const hi = [-Infinity, -Infinity, -Infinity];
          let largest = 0;
          for (const i of instances) {
            for (let a = 0; a < 3; ++a) {
              lo[a] = Math.min(lo[a], spheres[i * 4 + a]);
              hi[a] = Math.max(hi[a], spheres[i * 4 + a]);
            }
            largest = Math.max(largest, spheres[i * 4 + 3]);
          }
          for (let a = 0; a < 3; ++a) expect(gpuSpheres[d * 4 + a]).toBeCloseTo((lo[a] + hi[a]) / 2, 3);
          const expected = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) / 2 + largest;
          expect(r).toBeCloseTo(expected, 3);
        }
        for (const b of [data, segmentData, instanceCells, cellSpheres]) b.ref.value.destroy();
      }
      builder.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
