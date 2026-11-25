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

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
END
`;

/** 1cbs, or `small`: one residue */
async function addSpacefill(plugin: HeadlessPluginContext, small = false) {
  const data = await plugin.builders.data.rawData(
    small ? { data: Pdb, label: 'test' } : { data: readExample('1cbs_updated.cif'), label: '1cbs' },
  );
  const trajectory = await plugin.builders.structure.parseTrajectory(data, small ? 'pdb' : 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  // the structure's spheres as one plain object, which the raster takes with 'on' (a merged one only culled on the GPU)
  await plugin.builders.structure.representation.addRepresentation(structure, {
    type: 'spacefill',
    typeParams: { visuals: ['structure-element-sphere'] } as never,
  });
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
// pass's requests), so the first frame after the commit draws the raster's spheres instead of the impostors. Each case
// on new devices (the first frame's compiles are what it tests); medium: 1cbs, a second or more a case
describe.skipIf(!WebGPU.module)('sphere raster prepared on WebGPU', () => {
  // quick: one residue, prepared
  it('while representations build: drawn in the first frame, one residue', async () => {
    const { plugin, frame } = await createPlugin('on');
    try {
      const canvas = plugin.canvas3d!;
      if (!canvas.graphics.capabilities.atomicMinMax64) return;
      frame();
      await canvas.graphics.resources.pending.whenReady();
      await addSpacefill(plugin, true);
      // the resolve drawn: its pipeline and the kernels were ready (else the impostors draw that frame)
      const { compiled, rasterDrawn } = await firstFrame(plugin, frame);
      expect(compiled).toEqual([]);
      expect(rasterDrawn).toBe(true);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  for (const prepare of [true, false]) {
    it(
      `while representations build: ${prepare ? 'drawn' : 'without it, not drawn'} in the first frame`,
      { tags: ['medium'] },
      async () => {
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
      },
      120_000,
    );
  }

  it(
    'by a recorded precompile list (its stand-ins), before the objects exist',
    { tags: ['medium'] },
    async () => {
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
    },
    120_000,
  );
});
