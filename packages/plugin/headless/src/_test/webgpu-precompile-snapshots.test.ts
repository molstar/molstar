/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginState } from '@molstar/plugin/state';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { OverpaintStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { WebGPU } from './webgpu-module.js';

const readExample = (name: string) =>
  readFileSync(new URL(`../../../../../data/examples/${name}`, import.meta.url), 'utf8');

const createPlugin = async () => {
  const plugin = await HeadlessPluginContext.create(
    { webgpu: WebGPU.module! },
    DefaultPluginSpec(),
    { width: 32, height: 32 },
    { backend: 'webgpu' },
  );
  await plugin.init();
  return plugin;
};

const Pdb = `ATOM      1  N   ALA A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       0.551   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1      -0.200   2.390   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       0.540  -0.760  -1.210  1.00  0.00           C
END
`;

/** A state of two snapshots of one residue: the structure, then its spacefill as well */
async function createSmallState() {
  // no canvas (nor camera): nothing compiled for the state
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  try {
    const data = await plugin.builders.data.rawData({ data: Pdb, label: 'test' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    const { snapshot } = plugin.managers;
    snapshot.add(PluginStateSnapshotManager.Entry(plugin.state.getSnapshot({ camera: false }), {}));
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
    snapshot.add(PluginStateSnapshotManager.Entry(plugin.state.getSnapshot({ camera: false }), {}));
    // as written to and read from a file
    return JSON.parse(JSON.stringify(await snapshot.getStateSnapshot())) as PluginStateSnapshotManager.StateSnapshot;
  } finally {
    plugin.dispose();
  }
}

/** A state of two snapshots: a cartoon, then surfaces, transparent ball-and-stick and overpainted spacefill as well */
async function createState() {
  const plugin = await createPlugin();
  try {
    const data = await plugin.builders.data.rawData({ data: readExample('1cbs_updated.cif'), label: '1cbs' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model);
    const { representation } = plugin.builders.structure;
    const { snapshot } = plugin.managers;
    await representation.addRepresentation(structure, { type: 'cartoon' });
    snapshot.add(PluginStateSnapshotManager.Entry(plugin.state.getSnapshot(), {}));

    await representation.addRepresentation(structure, { type: 'ball-and-stick', typeParams: { alpha: 0.5 } });
    await representation.addRepresentation(structure, { type: 'molecular-surface', color: 'hydrophobicity' });
    await representation.addRepresentation(structure, { type: 'gaussian-surface', typeParams: { xrayShaded: true } });
    const spacefill = await representation.addRepresentation(structure, { type: 'spacefill' });
    await plugin.build().to(spacefill).apply(OverpaintStructureRepresentation3DFromScript).commit();
    snapshot.add(PluginStateSnapshotManager.Entry(plugin.state.getSnapshot(), {}));
    // as written to and read from a file
    return JSON.parse(JSON.stringify(await snapshot.getStateSnapshot())) as PluginStateSnapshotManager.StateSnapshot;
  } finally {
    plugin.dispose();
  }
}

/** The render pipelines applying `snapshot` (and preparing its frames and first interaction) compiles */
async function getCompiledPipelines(plugin: HeadlessPluginContext, snapshot: PluginState.Snapshot) {
  const pending = plugin.canvas3d!.graphics.resources.pending;
  await pending.whenReady();
  const before = pending.trackedCount;
  await plugin.state.setSnapshot(snapshot);
  plugin.canvas3d!.commit(true);
  await plugin.canvas3d!.prepare({ interaction: true });
  await pending.whenReady();
  // render pipelines only: the GPU surfaces' compute kernels and readbacks are their own (prepared as they compute)
  return pending.trackedSince(before).filter((label) => label.startsWith('pipeline '));
}

// Quick: one residue, its state made without a canvas, one device; medium: 1cbs in five representations, three devices
// compiling their pipelines
describe.skipIf(!WebGPU.module)('states loaded on WebGPU', () => {
  it('prepare what their snapshots draw, so that applying the next one compiles nothing: one residue', async () => {
    const state = await createSmallState();
    expect(state.entries.length).toBe(2);
    const plugin = await createPlugin();
    try {
      // the impostors alone (the raster's kernels are their own)
      plugin.canvas3d!.setProps({ renderer: { sphereRaster: 'off' } });
      // loading applies the first snapshot and prepares both in the background
      await plugin.managers.snapshot.setStateSnapshot(state);
      expect(await getCompiledPipelines(plugin, state.entries[1].snapshot)).toEqual([]);
    } finally {
      plugin.dispose();
    }
  }, 60_000);

  it(
    'prepare what their snapshots draw, so that applying the next one compiles nothing',
    { tags: ['medium'] },
    async () => {
      const state = await createState();
      expect(state.entries.length).toBe(2);
      const [, second] = state.entries;

      const a = await createPlugin();
      const b = await createPlugin();
      try {
        // loading applies the first snapshot and prepares both in the background
        await a.managers.snapshot.setStateSnapshot(state);
        expect(await getCompiledPipelines(a, second.snapshot)).toEqual([]);

        // control: the first snapshot applied without loading the state, the second compiles
        await b.state.setSnapshot(state.entries[0].snapshot);
        expect((await getCompiledPipelines(b, second.snapshot)).length).toBeGreaterThan(0);
      } finally {
        a.dispose();
        b.dispose();
      }
    },
    120_000,
  );
});
