/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { asWebGPUTexture } from '@molstar/graphics/gl/webgpu/texture';
import { WebGPU } from './webgpu-module.js';

// pass targets take memory only once drawn into: `define` records the size, the first render pass, copy or read
// allocates (unwritten textures are sampled as a placeholder)
describe.skipIf(!WebGPU.module)('lazy texture storage on WebGPU', () => {
  it('allocates render targets when first drawn into and again after a resize', async () => {
    const gpu = WebGPU.module!.create([]);
    const device = await requestWebGPUDevice('default', gpu);
    if (!device) throw new Error('no WebGPU device');
    const graphics = createWebGPUContext(device, { width: 16, height: 16 });
    try {
      const { resources } = graphics;
      const target = graphics.createRenderTarget(64, 32, 'depth-stencil');
      const texture = asWebGPUTexture(target.texture);
      // defined, not allocated
      expect(texture.getWidth()).toBe(64);
      expect(texture.hasStorage).toBe(false);
      expect(texture.getView()).toBeUndefined();
      expect(target.getByteCount()).toBe(0);
      expect(resources.getByteCounts().texture).toBe(0);
      expect(resources.getByteCounts().renderbuffer).toBe(0);

      graphics.beginRenderPass(target, { clear: { color: [0, 0, 0, 0], depth: 0 } });
      graphics.flush();
      expect(texture.hasStorage).toBe(true);
      expect(texture.getView()).toBeDefined();
      expect(target.getByteCount()).toBe(64 * 32 * 4 + target.depthRenderbuffer!.getByteCount());
      expect(target.depthRenderbuffer!.getByteCount()).toBeGreaterThan(0);

      // a new size releases the storage until drawn into again
      target.setSize(128, 64);
      expect(texture.hasStorage).toBe(false);
      expect(target.getByteCount()).toBe(0);
      graphics.beginRenderPass(target, { clear: { color: [0, 0, 0, 0], depth: 0 } });
      graphics.flush();
      expect(target.getByteCount()).toBeGreaterThan(128 * 64 * 4);

      // read before written: allocated, zeros
      const unwritten = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      unwritten.define(4, 4);
      expect(unwritten.getByteCount()).toBe(0);
      const data = await graphics.readTextureAsync(unwritten, 0, 0, 4, 4);
      expect(Array.from(data).every((v) => v === 0)).toBe(true);
      expect(unwritten.getByteCount()).toBe(4 * 4 * 4);

      target.destroy();
      unwritten.destroy();
      expect(graphics.takeErrors?.() ?? []).toEqual([]);
    } finally {
      graphics.destroy();
    }
    expect(gpu).toBeDefined();
  });
});
