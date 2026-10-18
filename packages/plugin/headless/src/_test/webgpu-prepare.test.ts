/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createCopyRenderable } from '@molstar/graphics/gl/compute/util';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { createWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { WebGPU } from './webgpu-module.js';

describe.skipIf(!WebGPU.module)('preparePipelines on WebGPU', () => {
  it('compiles the requested pipelines, so that draws with them start no work', async () => {
    const gpu = WebGPU.module!.create([]);
    const device = await requestWebGPUDevice('default', gpu);
    if (!device) throw new Error('no WebGPU device');
    const graphics = createWebGPUContext(device, { width: 16, height: 16 });
    try {
      const { pending } = graphics.resources;
      const target = graphics.createRenderTarget(16, 16);
      const copy = createCopyRenderable(graphics, graphics.createRenderTarget(16, 16).texture);
      const draw = (state: RenderState) => {
        graphics.beginRenderPass(target);
        graphics.setRenderState(state);
        copy.render();
        graphics.flush();
      };

      const descriptor = graphics.describeTarget(target);
      expect(descriptor).toEqual(graphics.describeTarget(target.framebuffer));
      expect(descriptor.key).not.toBe(graphics.describeTarget('drawing-buffer').key);

      await graphics.resources.preparePipelines([
        { item: copy.renderItem, variant: 'compute', state: RenderState.Default, target: descriptor },
      ]);

      // prepared: the draw finds its pipeline
      let tracked = pending.trackedCount;
      draw(RenderState.Default);
      expect(pending.trackedCount).toBe(tracked);

      // not prepared (another render state): the draw starts compiling its pipeline
      tracked = pending.trackedCount;
      draw(RenderState.Additive);
      expect(pending.trackedCount).toBe(tracked + 1);
      await pending.whenReady();
    } finally {
      graphics.destroy();
    }
    // keeps `gpu` alive until the context is gone (collecting it while its device is in use crashes Dawn's bindings)
    expect(gpu).toBeDefined();
  });
});
