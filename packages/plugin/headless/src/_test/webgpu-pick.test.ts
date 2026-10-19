/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Vec2 } from '@molstar/core/math/linear-algebra';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
END
`;

// a new canvas (as after a restore or a backend switch) picks before anything is drawn: the pick pipelines compile in
// the background and their draws are skipped, the pick waits for them instead of resolving to nothing
describe.skipIf(!WebGPU.module)('picking on WebGPU', () => {
  it('finds the structure with the first pick after loading', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 64, height: 64 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      const canvas = plugin.canvas3d!;
      // the scene committed, nothing drawn yet: no pipelines compiled
      canvas.commit(true);
      const pixelRatio = canvas.graphics.pixelRatio;
      const skippedDraws = canvas.graphics.stats.skippedDraws;

      const pick = await canvas.asyncIdentify(Vec2.create(32 / pixelRatio, 32 / pixelRatio))?.wait();
      // the case under test: the first pick render skipped its draws
      expect(canvas.graphics.stats.skippedDraws).toBeGreaterThan(skippedDraws);
      expect(pick?.id.objectId).toBeGreaterThanOrEqual(0);
    } finally {
      plugin.dispose();
    }
  });
});
