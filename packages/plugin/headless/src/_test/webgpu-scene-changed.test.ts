/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { Scene } from '@molstar/graphics/gl/scene';
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

describe.skipIf(!WebGPU.module)('Scene changed objects on WebGPU', () => {
  it('an object updated by `update(objects)` is not updated again by `updateChanged`', async () => {
    const plugin = await HeadlessPluginContext.create(
      { webgpu: WebGPU.module! },
      DefaultPluginSpec(),
      { width: 32, height: 32 },
      { backend: 'webgpu' },
    );
    try {
      await plugin.init();
      const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'cartoon' });
      await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
      plugin.canvas3d!.commit(true);

      const scene = Scene.create(plugin.canvas3d!.graphics);
      const objects = plugin.canvas3d!.getRenderObjects();
      for (const ro of objects) scene.add(ro);
      scene.commit();
      scene.updateChanged(false);
      const [a, b] = objects;
      const updateA = vi.spyOn(scene.get(a)!, 'update');
      const updateB = vi.spyOn(scene.get(b)!, 'update');

      // both changed (a recolor), `a` updated on its own (as `Canvas3D.update(repr)`), then a hover's update
      ValueCell.update(a.values.uMetalness, 0.3);
      ValueCell.update(b.values.uMetalness, 0.3);
      scene.update([a], false);
      expect(updateA).toHaveBeenCalledTimes(1);
      scene.updateChanged(false);
      expect(updateA).toHaveBeenCalledTimes(1);
      expect(updateB).toHaveBeenCalledTimes(1);

      // changed again after its own update: taken again
      ValueCell.update(a.values.uMetalness, 0.4);
      scene.updateChanged(false);
      expect(updateA).toHaveBeenCalledTimes(2);
    } finally {
      plugin.dispose();
    }
  });
});
