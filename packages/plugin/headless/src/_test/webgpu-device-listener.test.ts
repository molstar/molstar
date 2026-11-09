/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { createWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// A destroyed WebGPU context leaves nothing of itself on its device: Dawn's Node.js bindings keep the device by a
// handle of their own, also destroyed, and the error listener's scope was the whole context (every context of a test
// run held to its end: a software adapter's memory, in RAM, ran out)
describe.skipIf(!WebGPU.module)('a destroyed WebGPU context', () => {
  it('removes its error listener from its device', async () => {
    // the module's globals (usage flags) set up as for the plugin
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      const device = (await requestWebGPUDevice('default', WebGPU.module!.create([])))!;
      const added = vi.spyOn(device, 'addEventListener');
      const removed = vi.spyOn(device, 'removeEventListener');
      const ctx = createWebGPUContext(device, { width: 16, height: 16, ownsDevice: true });
      const listeners = added.mock.calls.filter(([type]) => type === 'uncapturederror').map(([, l]) => l);
      expect(listeners).toHaveLength(1);
      ctx.destroy();
      expect(removed.mock.calls.filter(([type]) => type === 'uncapturederror').map(([, l]) => l)).toEqual(listeners);
    } finally {
      plugin.dispose();
    }
  });
});
