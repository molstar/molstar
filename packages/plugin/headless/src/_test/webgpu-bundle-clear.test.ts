/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Size = 8;

// clears between the draws of a recorded pass (e.g. the stencil of solid interior caps) are draws of the render bundle,
// replayed as cleared; where the bundle can't clear as the pass would, the recording fails (drawn directly)
describe.skipIf(!WebGPU.module)('clears in render bundles on WebGPU', () => {
  it('replay as cleared, or fail the recording', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: Size, height: Size },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const g = plugin.canvas3d!.graphics as WebGPUContext;
      const target = g.createRenderTarget(Size, Size, 'depth-stencil', 'uint8', 'nearest');
      const begin = () => {
        g.beginRenderPass(target, { label: 'test', clear: { color: [0, 0, 1, 1], depth: 0, stencil: 0 } });
        g.setViewport(0, 0, Size, Size);
        g.setScissor(0, 0, Size, Size);
      };
      const pixel = async () => Array.from((await g.readRenderTarget(target, 0, 0, 1, 1)).slice(0, 4));

      // a color clear after the pass began: recorded, replayed in a pass cleared otherwise
      begin();
      const bundle = g.recordRenderBundle(() => g.clear({ color: [1, 0, 0, 1], stencil: 0 }));
      expect(bundle).toBeDefined();
      g.executeRenderBundles([bundle!]);
      expect(await pixel()).toEqual([255, 0, 0, 255]);
      begin();
      expect(await pixel()).toEqual([0, 0, 255, 255]);
      begin();
      g.executeRenderBundles([bundle!]);
      expect(await pixel()).toEqual([255, 0, 0, 255]);

      // a stencil value other than 0 needs the reference (pass state): not recorded
      begin();
      expect(g.recordRenderBundle(() => g.clear({ stencil: 1 }))).toBeUndefined();
      // a scissor beyond the viewport: the bundle's clear covers the viewport only
      begin();
      g.setViewport(0, 0, Size / 2, Size / 2);
      expect(g.recordRenderBundle(() => g.clear({ color: [1, 0, 0, 1] }))).toBeUndefined();
      target.destroy();
    } finally {
      plugin.dispose();
    }
  });
});
