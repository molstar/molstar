/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// right after a real loss there may be no device yet (the GPU process restarting): `recreate` asks again
describe.skipIf(!WebGPU.module)('recreating a lost WebGPU context', () => {
  it('asks for a device again while none is given', async () => {
    // the module's globals (usage flags) set up as for the plugin
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    const gpu = WebGPU.module!.create([]);
    const device = await requestWebGPUDevice('default', gpu);
    let requests = 0;
    const ctx = createWebGPUContext(device!, {
      width: 16,
      height: 16,
      requestDevice: async () => (++requests < 3 ? undefined : requestWebGPUDevice('default', gpu)),
      ownsDevice: true,
    });
    try {
      const lost = new Promise<void>((resolve) => ctx.contextLost.subscribe(() => resolve()));
      expect(ctx.simulateContextLoss()).toBe(true);
      await lost;
      const next = await ctx.recreate();
      expect(requests).toBe(3);
      expect(next).toBeDefined();
      expect(next!.isContextLost).toBe(false);
      next!.destroy();
    } finally {
      ctx.destroy({ lost: true });
      plugin.dispose();
    }
  });
});
