/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { ReleaseDelayMs } from '@molstar/graphics/canvas3d/passes/pass-targets';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const Pdb = `${Array.from({ length: 20 }, (_, i) => {
  const f = (v: number) => v.toFixed(3).padStart(8);
  const serial = String(i + 1).padStart(5);
  const residue = String(i + 1).padStart(4);
  return `ATOM  ${serial}  CA  ALA A${residue}    ${f(Math.cos(i) * 8)}${f(Math.sin(i) * 8)}${f(i * 1.5)}  1.00  0.00           C`;
}).join('\n')}\nEND\n`;

const wait = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

// the headless canvas context takes its props as the browser's: transparency and pick scale reach the passes
describe.skipIf(!WebGPU.module)('headless canvas context props on WebGPU', () => {
  it('switch the transparency mode, holding frames until it is prepared, and the pick scale', async () => {
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
      await plugin.builders.structure.representation.addRepresentation(structure, {
        type: 'spacefill',
        typeParams: { alpha: 0.5 },
      });
      const canvas = plugin.canvas3d!;
      const context = plugin.renderer.canvas3dContext!;
      const { passes } = context;
      const frames = async () => {
        for (let i = 0; i < 2; ++i) {
          canvas.commit(true);
          canvas.requestDraw();
          canvas.tick(now(), { isSynchronous: true });
          await canvas.graphics.resources.pending.whenReady();
        }
      };

      await frames();
      expect(passes.draw.transparency).toBe('wboit');
      expect(passes.draw.wboit.getByteCount()).toBeGreaterThan(0);

      // the frames wait for the pipelines of the new mode, so that none draws with its draws skipped (an empty frame)
      const skippedDraws = canvas.graphics.stats.skippedDraws;
      plugin.setCanvas3DContextProps({ transparency: 'dpoit' });
      expect(context.props.transparency).toBe('dpoit');
      expect(passes.draw.transparency).toBe('dpoit');
      await frames();
      expect(canvas.graphics.stats.skippedDraws).toBe(skippedDraws);
      expect(passes.draw.dpoit.getByteCount()).toBeGreaterThan(0);
      // the targets of the mode no longer drawn are released
      await wait(ReleaseDelayMs + 200);
      expect(passes.draw.wboit.getByteCount()).toBe(0);
      // and an image renders with it
      const image = await plugin.getImageRaw({ width: 64, height: 64 });
      expect(image.width).toBe(64);

      const pickRatio = passes.pick.pickRatio;
      plugin.setCanvas3DContextProps({ pickScale: context.props.pickScale / 2 });
      expect(passes.pick.pickRatio).toBe(pickRatio / 2);
    } finally {
      plugin.dispose();
    }
  });
});
