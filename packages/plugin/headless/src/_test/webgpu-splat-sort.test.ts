/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createSplatSorter, SplatKeyCount, SplatSortKernels } from '@molstar/graphics/gl/compute/splat-sort';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

/** `count` centers (`vec4` each) of a fixed sequence (a linear congruential generator), within ±50 */
function centers(count: number) {
  const array = new Float32Array(count * 4);
  let x = 2024;
  for (let i = 0; i < count * 4; i++) {
    x = (Math.imul(x, 1103515245) + 12345) >>> 0;
    array[i] = (i & 3) === 3 ? 0 : (x / 2 ** 32) * 100 - 50;
  }
  return array;
}

describe.skipIf(!WebGPU.module)('gaussian splats sorted on the GPU', () => {
  // the order of two instances, each its block: every splat once, back to front by depth within a key's step
  it('sorts each instance in its block back to front, every splat once', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const compute = plugin.canvas3d!.graphics.compute!;
      await compute.prepareKernels(SplatSortKernels);
      const count = 200_000;
      const data = centers(count);
      const positions = compute.createBuffer(count * 16, 'centers');
      positions.write(data);
      const order = compute.createBuffer(count * 2 * 4, 'order');
      const sorter = createSplatSorter(compute);
      // the depth range: ±50 along a unit row, 87 along (1, 1, 1) / √3
      const rows = [
        [0, 0, 1],
        [1 / Math.sqrt(3), -1 / Math.sqrt(3), 1 / Math.sqrt(3)],
      ];
      const range = 50 * Math.sqrt(3);
      sorter.sort(
        positions,
        order,
        rows.map(([a, b, c], k) => ({ offset: k * count, count, a, b, c, depthMin: -range, depthMax: range })),
      );
      const words = new Uint32Array(await order.readAsync());
      const step = (2 * range) / (SplatKeyCount - 1);
      rows.forEach(([a, b, c], k) => {
        const block = words.subarray(k * count, (k + 1) * count);
        const seen = new Uint8Array(count);
        let previous = -Infinity;
        let outOfOrder = 0;
        for (const i of block) {
          seen[i]++;
          const z = a * data[i * 4] + b * data[i * 4 + 1] + c * data[i * 4 + 2];
          // within a key: in any order; a pair just over a step apart may share a key (the GPU's depths are f32, the
          // key's boundary rounded), and the atomics order a key's splats differently from run to run
          if (z < previous - step * 1.01) outOfOrder++;
          previous = Math.max(previous, z);
        }
        expect(seen.every((v) => v === 1)).toBe(true);
        expect(outOfOrder).toBe(0);
      });
      sorter.destroy();
      positions.destroy();
      order.destroy();
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
