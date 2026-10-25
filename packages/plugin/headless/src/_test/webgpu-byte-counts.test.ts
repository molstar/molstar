/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
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

describe.skipIf(!WebGPU.module)('resource counts on WebGPU', () => {
  it('count the live resources and their bytes, released with the scene', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const { resources, stats } = plugin.canvas3d!.graphics;
      const counts = () => ({
        ...stats.resourceCounts,
        drawCount: stats.drawCount,
        instanceCount: stats.instanceCount,
        instancedDrawCount: stats.instancedDrawCount,
      });
      const load = async () => {
        const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
        const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
        await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
        await plugin.getImageRaw({ width: 64, height: 64 });
        return resources.getByteCounts();
      };
      const clear = async () => {
        await plugin.clear();
        plugin.canvas3d!.commit(true);
        await plugin.getImageRaw({ width: 64, height: 64 });
        return resources.getByteCounts();
      };

      const drawn = await load();
      // data textures of the representations and the targets of the passes
      expect(drawn.texture).toBeGreaterThan(0);
      // cartoon meshes are indexed
      expect(drawn.attribute).toBeGreaterThan(0);
      expect(drawn.elements).toBeGreaterThan(0);
      expect(drawn.uniform).toBeGreaterThan(0);
      // no read in flight
      expect(drawn.pixelPack).toBe(0);
      const loaded = counts();
      expect(loaded.attribute).toBeGreaterThan(0);
      expect(loaded.elements).toBeGreaterThan(0);
      expect(loaded.uniform).toBeGreaterThan(0);
      expect(loaded.shader).toBeGreaterThan(0);
      expect(loaded.pixelPack).toBe(0);
      expect(loaded.instancedDrawCount).toBeGreaterThan(0);

      const cleared = await clear();
      // what remains is the helpers' geometry (e.g. the camera helper), far less than the structure's
      expect(cleared.attribute).toBeLessThan(drawn.attribute);
      // what remains of elements is the context's shared quad indices (`quad-indices.ts`), kept for the next quads
      expect(cleared.elements).toBeLessThan(drawn.elements);
      // texture bytes include the pass targets, allocated by what frames draw (an empty scene draws other passes), so
      // the totals of two scenes don't compare; the structure's textures are released by the check below
      // the structure's resources are released: after helpers and uniform arenas settled (a first cycle), loading and
      // clearing again leaves the same counts
      await load();
      const settled = await clear();
      const settledCounts = counts();
      await load();
      expect(await clear()).toEqual(settled);
      // and the same number of live resources (counted until destroyed) and draws
      expect(counts()).toEqual(settledCounts);
    } finally {
      plugin.dispose();
    }
  });
});
