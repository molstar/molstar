/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `${Array.from({ length: 16 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  const [x, y, z] = [(i % 4) * 2.5, Math.floor(i / 4) * 2.5, 0];
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(x)}${f(y)}${f(z)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

// an update of a representation whose new defines need programs not compiled yet (here another color type) holds the
// frames until they are ready, as a props change does: the last frame stays shown, instead of frames without the
// representation (its draws skipped while the pipelines compile)
describe.skipIf(!WebGPU.module)('frames after a representation update on WebGPU', () => {
  it('wait for the pipelines the update needs instead of skipping its draws', async () => {
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
      const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        color: 'element-symbol',
      });
      const canvas = plugin.canvas3d!;
      const { stats } = canvas.graphics;
      const { pending } = canvas.graphics.resources;
      let time = now();
      let drawn = 0;
      const sub = canvas.didDraw.subscribe(() => (drawn += 1));
      const frame = () => {
        canvas.requestDraw();
        canvas.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
      };
      frame();
      await pending.whenReady();
      frame();
      await pending.whenReady();

      // another color type: programs (and pipelines) the representation has not drawn with yet
      await plugin
        .build()
        .to(repr)
        .update((old) => ({ ...old, colorTheme: { name: 'uniform', params: {} } }))
        .commit();
      const skipped = stats.skippedDraws;
      const drawnBefore = drawn;
      frame();
      // the pipelines compile in the background meanwhile: no frame drawn, none skipping draws
      expect(pending.trackedCount).toBeGreaterThan(0);
      expect(stats.skippedDraws).toBe(skipped);
      expect(drawn).toBe(drawnBefore);
      await pending.whenReady();
      frame();
      expect(stats.skippedDraws).toBe(skipped);
      expect(drawn).toBeGreaterThan(drawnBefore);
      sub.unsubscribe();
    } finally {
      plugin.dispose();
    }
  }, 60_000);
});
