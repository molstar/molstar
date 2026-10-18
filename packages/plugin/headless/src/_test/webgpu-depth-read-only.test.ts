/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { createCopyRenderable } from '@molstar/graphics/gl/compute/util';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 4;

describe.skipIf(!WebGPU.module)('read-only depth passes on WebGPU', () => {
  it('bind the depth attachment as a texture, which a writable one can not', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const graphics = plugin.canvas3d!.graphics as WebGPUContext;
      const { resources } = graphics;
      const color = resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
      color.define(Size, Size);
      const depth = resources.texture('image-depth', 'depth-stencil', 'float-stencil', 'nearest');
      depth.define(Size, Size);
      const framebuffer = resources.framebuffer();
      color.attachFramebuffer(framebuffer, 'color0');
      depth.attachFramebuffer(framebuffer, 'depth');
      // the depth in the red channel of the color target
      const copy = createCopyRenderable(graphics, depth);
      const pass = (depthReadOnly: boolean) => {
        graphics.beginRenderPass(framebuffer, { clear: { color: [0, 0, 0, 0] } });
        graphics.endRenderPass();
        graphics.beginRenderPass(framebuffer, { depthReadOnly });
        graphics.setViewport(0, 0, Size, Size);
        copy.render();
        graphics.endRenderPass();
      };
      // the first draw of a pipeline starts compiling it in the background (skipped), the second draws
      const draw = async (depthReadOnly: boolean) => {
        pass(depthReadOnly);
        await resources.pending.whenReady();
        pass(depthReadOnly);
        await graphics.waitForGpuCommandsComplete();
      };

      graphics.beginRenderPass(framebuffer, { clear: { depth: 0.25, stencil: 0 } });
      graphics.endRenderPass();
      await draw(true);
      const pixels = await graphics.readTextureData(color);
      expect(Math.abs(pixels.array[0] - 0.25 * 255)).toBeLessThanOrEqual(1);
      expect(graphics.takeErrors()).toEqual([]);

      // as an attachment that is written, the texture can not also be read in the pass
      await draw(false);
      expect(graphics.takeErrors().length).toBeGreaterThan(0);

      graphics.beginRenderPass(framebuffer, { depthReadOnly: true });
      expect(() => graphics.clear({ depth: 0 })).toThrow(/read-only/);
      graphics.endRenderPass();
    } finally {
      plugin.dispose();
    }
  });
});
