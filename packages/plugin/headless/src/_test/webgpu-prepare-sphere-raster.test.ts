/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { now } from '@molstar/core/util/now';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) =>
  readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');

async function createPlugin(sphereRaster: 'on' | 'off') {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: 64, height: 64 },
    { backend: 'webgpu' },
  );
  await plugin.init();
  plugin.canvas3d!.setProps({
    renderer: { sphereRaster, sphereRasterRadius: 32 },
    camera: { mode: 'perspective' } as never,
    postprocessing: { occlusion: { name: 'off', params: {} }, outline: { name: 'off', params: {} } },
  });
  let time = now();
  const frame = () => {
    plugin.canvas3d!.requestDraw();
    plugin.canvas3d!.tick((time += 1000) as now.Timestamp, { isSynchronous: true });
  };
  return { plugin, frame };
}

async function addSpacefill(plugin: HeadlessPluginContext) {
  const data = await plugin.builders.data.rawData({ data: readExample('1cbs_updated.cif'), label: '1cbs' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
}

/** The first frame after the objects were added: the pipelines it compiled, whether the raster drew */
async function firstFrame(plugin: HeadlessPluginContext, frame: () => void) {
  const { graphics } = plugin.canvas3d!;
  const { pending } = graphics.resources;
  // what was prepared is ready by now
  await pending.whenReady();
  const resolves = graphics.stats.sphereRasterResolves;
  const before = pending.trackedCount;
  frame();
  const compiled = pending.trackedSince(before).filter((label) => label.startsWith('pipeline '));
  const rasterDrawn = graphics.stats.sphereRasterResolves > resolves;
  await pending.whenReady();
  expect(graphics.takeErrors?.() ?? []).toEqual([]);
  return { compiled, rasterDrawn };
}

// with the sphere raster on, what is prepared covers the raster: its kernels and the resolve's pipeline (the opaque
// pass's requests), so the first frame after the commit draws the raster's spheres instead of the impostors
describe.skipIf(!WebGPU.module)('sphere raster prepared on WebGPU', () => {
  for (const prepare of [true, false]) {
    it(`while representations build: ${prepare ? 'drawn' : 'without it, not drawn'} in the first frame`, async () => {
      const { plugin, frame } = await createPlugin('on');
      try {
        const canvas = plugin.canvas3d!;
        if (!canvas.graphics.capabilities.atomicMinMax64) return;
        if (!prepare) canvas.prepareObjects = () => {};
        frame();
        await canvas.graphics.resources.pending.whenReady();
        await addSpacefill(plugin);
        const { compiled, rasterDrawn } = await firstFrame(plugin, frame);
        // the resolve drawn: its pipeline and the kernels were ready (else the impostors draw that frame)
        if (prepare) expect(compiled).toEqual([]);
        expect(rasterDrawn).toBe(prepare);
      } finally {
        plugin.dispose();
      }
    }, 120_000);
  }

  it('by a recorded precompile list (its stand-ins), before the objects exist', async () => {
    const a = await createPlugin('off');
    const b = await createPlugin('on');
    try {
      if (!b.plugin.canvas3d!.graphics.capabilities.atomicMinMax64) return;
      const canvasA = a.plugin.canvas3d!;
      canvasA.resetPrecompileList();
      await addSpacefill(a.plugin);
      // recorded at the commit
      canvasA.commit(true);
      const list = JSON.parse(JSON.stringify(canvasA.getPrecompileList()));
      expect(list.objects.some((o: { type: string }) => o.type === 'spheres')).toBe(true);

      const canvasB = b.plugin.canvas3d!;
      // the list under test only
      canvasB.prepareObjects = () => {};
      b.frame();
      await canvasB.graphics.resources.pending.whenReady();
      await canvasB.precompile(list);
      await addSpacefill(b.plugin);
      // the resolve drawn: its pipeline and the kernels were ready (else the impostors draw that frame)
      const { compiled, rasterDrawn } = await firstFrame(b.plugin, b.frame);
      expect(rasterDrawn).toBe(true);
      expect(compiled).toEqual([]);
    } finally {
      a.plugin.dispose();
      b.plugin.dispose();
    }
  }, 120_000);
});
