/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { HeadlessPluginContext } from '@molstar/plugin-headless/context';
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

describe.skipIf(!WebGPU.module)('states loaded on WebGPU', () => {
  it('prepare what their snapshots draw, so that applying the next one compiles nothing', async () => {
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
  }, 120_000);
});
