/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import type { WebGPUContext } from '@molstar/graphics/gl/webgpu/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
ATOM      6  N   GLY A   2       1.870   1.550   0.000  1.00  0.00           N
ATOM      7  CA  GLY A   2       2.520   2.860   0.000  1.00  0.00           C
ATOM      8  C   GLY A   2       4.030   2.720   0.000  1.00  0.00           C
ATOM      9  O   GLY A   2       4.570   1.610   0.000  1.00  0.00           O
END
`;

describe.skipIf(!WebGPU.module)('HeadlessPluginContext on WebGPU', () => {
  it('renders a structure', async () => {
    const webgpu = WebGPU.module!;

    const plugin = await HeadlessPluginContext.create(
      { webgpu },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      expect(plugin.canvas3d!.graphics.capabilities.backend).toBe('webgpu');

      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');

      const image = await plugin.getImageRaw({ width: 64, height: 64 });
      expect(image.data.length).toBe(64 * 64 * 4);
      // the structure covers part of the image, the rest is background
      const first = image.data.slice(0, 4).join();
      let differing = 0;
      for (let i = 0; i < image.data.length; i += 4) {
        if (image.data.slice(i, i + 4).join() !== first) differing += 1;
      }
      expect(differing).toBeGreaterThan(64);
      expect(differing).toBeLessThan(64 * 64);

      // device loss: a new context with a new device replaces the lost one, Canvas3D is built again on it, the image
      // is the same
      const graphics = plugin.canvas3d!.graphics as WebGPUContext;
      const lostDevice = graphics.device;
      const restored = new Promise<void>((resolve) =>
        plugin.canvas3d!.contextChanged.subscribe(({ reason }) => {
          if (reason === 'restore') resolve();
        }),
      );
      lostDevice.destroy();
      await restored;
      const next = plugin.canvas3d!.graphics as WebGPUContext;
      expect(next).not.toBe(graphics);
      expect(next.device).not.toBe(lostDevice);
      expect(next.isContextLost).toBe(false);
      const restoredImage = await plugin.getImageRaw({ width: 64, height: 64 });
      expect(restoredImage.data).toEqual(image.data);

      // release on purpose (e.g. a canvas out of view): the device is destroyed and the loss reported, but the context
      // stays lost until `restore`, which builds everything again on a new device, the image the same
      const context = plugin.renderer.canvas3dContext!;
      const lost = new Promise<void>((resolve) => context.contextLost!.subscribe(() => resolve()));
      let restores = 0;
      const sub = plugin.canvas3d!.contextChanged.subscribe(({ reason }) => {
        if (reason === 'restore') restores += 1;
      });
      expect(await context.restore!()).toBe(false);
      expect(context.release!()).toBe(true);
      expect(context.release!()).toBe(false);
      await lost;
      expect(context.isReleased).toBe(true);
      expect(context.graphics).toBe(next);
      expect(next.isContextLost).toBe(true);
      // not restored on its own
      await new Promise((resolve) => setTimeout(resolve, 100));
      expect(restores).toBe(0);
      expect(context.graphics).toBe(next);

      expect(await context.restore!()).toBe(true);
      sub.unsubscribe();
      expect(restores).toBe(1);
      expect(context.isReleased).toBe(false);
      const reacquired = plugin.canvas3d!.graphics as WebGPUContext;
      expect(reacquired).not.toBe(next);
      expect(reacquired.isContextLost).toBe(false);
      const reacquiredImage = await plugin.getImageRaw({ width: 64, height: 64 });
      expect(reacquiredImage.data).toEqual(image.data);
    } finally {
      plugin.dispose();
    }
  });
});
