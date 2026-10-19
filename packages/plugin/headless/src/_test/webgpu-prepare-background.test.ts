/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { BackgroundPass } from '@molstar/graphics/canvas3d/passes/background';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

// `Canvas3D.prepare` with props that may never be set (a snapshot's): the live background (its images) stays
describe.skipIf(!WebGPU.module)('preparing with other background props on WebGPU', () => {
  it('plans the background of the props without updating the live one', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 16, height: 16 },
      { backend: 'webgpu' },
    );
    const update = vi.spyOn(BackgroundPass.prototype, 'update');
    const prepare = vi.spyOn(BackgroundPass.prototype, 'prepare');
    try {
      await plugin.init();
      const canvas = plugin.canvas3d!;
      const { postprocessing } = canvas.props;
      update.mockClear();
      await canvas.prepare({
        props: {
          postprocessing: {
            ...postprocessing,
            enabled: true,
            background: {
              variant: {
                name: 'radialGradient',
                params: { centerColor: 0xffffff, edgeColor: 0x000000, ratio: 0.5, coverage: 'viewport' },
              },
            },
          },
        } as Parameters<typeof canvas.prepare>[0] extends { props?: infer P } ? P : never,
      });
      expect(update).not.toHaveBeenCalled();
      expect(prepare).toHaveBeenCalledWith(
        expect.anything(),
        expect.objectContaining({ variant: expect.objectContaining({ name: 'radialGradient' }) }),
      );
    } finally {
      update.mockRestore();
      prepare.mockRestore();
      plugin.dispose();
    }
  });
});
