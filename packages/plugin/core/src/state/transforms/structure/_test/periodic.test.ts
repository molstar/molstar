/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import type { Model, Structure } from '@molstar/model/model/structure';
import type { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { createPlaybackPlugin } from '../../../animation/_test/harness.js';

// The fixture: data/examples/trajectory/periodic/split-peptide.gro, a peptide split across a dodecahedral box in its
// first frame (between residues 2 and 3, a false polymer gap where it was read).
const gro = fs.readFileSync(
  path.resolve(import.meta.dirname, '../../../../../../../../data/examples/trajectory/periodic/split-peptide.gro'),
  'utf8',
);

async function load(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: gro, label: 'split-peptide.gro' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'gro');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { model, structure };
}

const setModel = (plugin: PluginContext, ref: string, params: { modelIndex: number; periodic?: 'whole' | 'off' }) =>
  plugin
    .build()
    .to(ref)
    .update(ModelFromTrajectory, () => params)
    .commit();

/** The CA–CA distance of residues 2 and 3, across the face in the first frame as read. */
function caca(model: Model) {
  const { x, y, z } = model.atomicConformation;
  // five atoms a residue: N CA C O CB
  return Math.hypot(x[6] - x[11], y[6] - y[11], z[6] - z[11]);
}

const gaps = (s: Structure) => s.units.reduce((n, u) => n + u.gapElements.length, 0);

describe('ModelFromTrajectory in a periodic box', () => {
  it('makes the molecules whole by default, and shows the file as read when asked', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model, structure } = await load(plugin);
    const cell = () => plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    const units = () => plugin.state.data.cells.get(structure.ref)!.obj!.data as Structure;

    expect(caca(cell())).toBeLessThan(4);
    expect(gaps(units())).toBe(0);

    await setModel(plugin, model.ref, { modelIndex: 0, periodic: 'off' });
    expect(caca(cell())).toBeGreaterThan(10);
    // the structure was remapped onto the frame as read: its false gap is back
    expect(gaps(units())).toBeGreaterThan(0);

    await setModel(plugin, model.ref, { modelIndex: 0, periodic: 'whole' });
    expect(caca(cell())).toBeLessThan(4);
    expect(gaps(units())).toBe(0);
  });

  it('keeps every frame whole, the same model for the same frame', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await load(plugin);
    const first = plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    for (let i = 1; i < 8; i++) {
      // as the animation updates it: the model index only
      await setModel(plugin, model.ref, { modelIndex: i });
      expect(caca(plugin.state.data.cells.get(model.ref)!.obj!.data as Model)).toBeLessThan(4);
    }
    await setModel(plugin, model.ref, { modelIndex: 0 });
    expect(plugin.state.data.cells.get(model.ref)!.obj!.data).toBe(first);
  });
});
