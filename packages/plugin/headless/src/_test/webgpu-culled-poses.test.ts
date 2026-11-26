/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { GridPositions } from '@molstar/core/math/geometry/grid-positions';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { ComputeKernel } from '@molstar/graphics/gl/compute';
import { createCulledPoseExpander } from '@molstar/graphics/gl/compute/culled-poses';
import { createU32RunWriter } from '@molstar/graphics/gl/compute/u32-runs';
import { getDispatchSize } from '@molstar/graphics/gl/compute/prefix-sum';
import { CulledData } from '@molstar/graphics/gl/culled-data';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createWebGPUPlugin } from './webgpu-plugin.js';
import { WebGPU } from './webgpu-module.js';

/** a particle list of `count` elements, a scale not above 0 among them (read as 1) */
function createList(count: number) {
  const positions = new Float32Array(count * 3);
  const rotations = new Float32Array(count * 4);
  const scales = new Float32Array(count);
  const q = Quat();
  for (let e = 0; e < count; ++e) {
    positions.set([e * 1.5 - 20, 40 - e, e * 0.25], e * 3);
    Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(1, e % 7, 2)), e * 0.9);
    Quat.toArray(q, rotations, e * 4);
    scales[e] = e % 13 === 3 ? 0 : 0.5 + (e % 10) * 0.1;
  }
  return { positions, rotations, scales };
}

describe.skipIf(!WebGPU.module)('poses of culled data expanded on the GPU', () => {
  it('writes the slots the CPU writes, for every kind of positions, rotations and scales, in batches; the runs too', async () => {
    const plugin = await createWebGPUPlugin(DefaultPluginSpec(), { width: 16, height: 16 });
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics;
      const compute = graphics.compute!;
      const list = createList(400);
      const snorm = Int16Array.from(list.rotations, (v) => Math.round(v * 32767));
      // a snorm16 view at an odd pair: its slices copied to be aligned
      const shifted = new Int16Array(snorm.length + 2);
      shifted.set(snorm, 2);
      const snormOdd = shifted.subarray(2);
      const grid = GridPositions.fromFloat(list.positions);
      const center = Vec3.create(0.5, -1, 2);
      // a similarity, folded into the poses: rotated, scaled uniformly, moved
      const similar = Mat4.mul(
        Mat4(),
        Mat4.fromTranslation(Mat4(), Vec3.create(3, -2, 1)),
        Mat4.mul(Mat4(), Mat4.fromRotation(Mat4(), 1.1, Vec3.create(0, 1, 0)), Mat4.fromUniformScaling(Mat4(), 2)),
      );
      const matrices = new Float32Array(5 * 16);
      for (let i = 0; i < 5; ++i) Mat4.toArray(Mat4.fromTranslation(Mat4(), Vec3.create(i, 2, -i)), matrices, i * 16);
      const indices = Array.from({ length: 60 }, (_, i) => (i * 37) % 400);
      const transforms = [
        InstanceTransforms.ofMatrices(matrices, 5),
        // a range in the grid's order, with the similarity
        InstanceTransforms.ofPoses({
          ...list,
          count: 50,
          first: 20,
          order: [...Array(50).keys()].reverse(),
          center,
          matrix: similar,
        }),
        // by indices, gathered
        InstanceTransforms.ofPoses({ ...list, count: 60, indices, center }),
        // on a grid, snorm16 (at an odd pair), a uniform scale
        InstanceTransforms.ofPoses({
          positions: grid,
          rotations: snormOdd,
          scale: 1.5,
          count: 120,
          first: 100,
          center,
        }),
        // no rotations, no scales; by indices on a grid
        InstanceTransforms.ofPoses({ positions: grid, count: 30, indices: indices.slice(0, 30), center }),
        InstanceTransforms.ofPoses({
          positions: list.positions,
          rotations: snorm,
          scales: list.scales,
          count: 150,
          first: 250,
          center,
        }),
      ];
      let first = 0;
      const parts = transforms.map((t) => {
        const part = { transforms: t, first };
        first += t.count;
        return part;
      });
      const n = first;
      const bases: number[] = [];
      for (const p of parts) bases.push(p.first);
      const segments = { count: parts.length, instanceBases: bases, instanceCounts: transforms.map((t) => t.count) };
      // serial but in the third segment (its words written by the CPU), the slot words of the second (in its own order)
      const aInstance = new Uint32Array(n).map((_, i) => i);
      aInstance.subarray(bases[2], bases[2] + transforms[2].count).reverse();
      const layout = CulledData.layout(parts, aInstance, segments, undefined, n, 64);
      expect(layout.poseSlots.length).toBe(5);
      expect(layout.runs.length).toBeGreaterThan(segments.count);

      const bytes = layout.vec4s * 16;
      const expected = new Uint8Array(bytes);
      layout.write((d, o) => expected.set(new Uint8Array(d.buffer, d.byteOffset, d.byteLength), o));
      const dispatchItems = (kernel: ComputeKernel, count: number) => {
        // rows of workgroups: the kernel linearizes the id
        const { x, y } = getDispatchSize(Math.ceil(count / 64), 3);
        kernel.dispatch(x, y);
      };
      // one batch; batches of about 40 ids, the parts split across them (each batch a submit)
      for (const batchBytes of [1 << 20, 40 * 5 * 4]) {
        const expander = createCulledPoseExpander(compute, batchBytes);
        const runWriter = createU32RunWriter(compute);
        const data = compute.createBuffer(bytes);
        // the CPU's without the poses and the runs first (a new buffer zeroed), the GPU's runs and poses after
        layout.write((d, o) => data.write(d, o), { poses: false, runs: false });
        runWriter.write(data, layout.runs, dispatchItems);
        expander.expand(layout, data, dispatchItems);
        const gpu = new Float32Array(await data.readAsync());
        const cpu = new Float32Array(expected.buffer);
        let differing = 0;
        for (let i = 0; i < cpu.length; ++i) {
          if (Math.abs(gpu[i] - cpu[i]) > 1e-4 * Math.max(1, Math.abs(cpu[i]))) differing += 1;
        }
        expect(differing).toBe(0);
        // the words before the slots (header, slot words, aInstance, segment indices) exact: as floats most are denormals
        const words = layout.slotFirst * 4;
        expect(Array.from(new Uint32Array(gpu.buffer, 0, words))).toEqual(
          Array.from(new Uint32Array(expected.buffer, 0, words)),
        );
        data.destroy();
        expander.destroy();
        runWriter.destroy();
      }
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
