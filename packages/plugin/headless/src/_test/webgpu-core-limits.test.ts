/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { InstanceCopyKernel, InstanceCopySchema } from '@molstar/graphics/gl/compute/cell-cull';
import { canCullOnGpu } from '@molstar/graphics/gl/renderable/merged';
import { createWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { WebGPU } from './webgpu-module.js';

// a device with the core limits (what every WebGPU device has): kernels binding more storage buffers than a stage may
// (8) are not supported, and so objects are culled on the CPU (the instance copy of GPU culling binds 10)
describe.skipIf(!WebGPU.module)('WebGPU core limits', () => {
  it('cull on the CPU where the culling kernels bind more storage buffers than a stage may', async () => {
    const gpu = WebGPU.module!.create([]);
    for (const limits of ['core', 'adapter'] as const) {
      const device = await requestWebGPUDevice('default', gpu, limits);
      if (!device) throw new Error('no WebGPU device');
      const graphics = createWebGPUContext(device, { width: 16, height: 16, ownsDevice: true });
      try {
        const fits = device.limits.maxStorageBuffersPerShaderStage >= 10;
        expect(graphics.compute!.isKernelSupported(InstanceCopyKernel, InstanceCopySchema)).toBe(fits);
        expect(canCullOnGpu(graphics)).toBe(fits && !!graphics.capabilities.indirectInstanceRanges);
        if (limits === 'core') expect(fits).toBe(false);
      } finally {
        graphics.destroy();
      }
    }
  });
});
