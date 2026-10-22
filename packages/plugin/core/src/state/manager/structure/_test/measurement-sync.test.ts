/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { StateSelection, StateTransform } from '@molstar/core/state';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { CreateGroup } from '@molstar/plugin/state/transforms/misc/group';
import { MeasurementGroupTag } from '../measurement.js';

const pdb = [
  'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  ',
  'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00  0.00           C  ',
].join('\n');

async function createPlugin() {
  const plugin = new PluginContext({ registry: DefaultPluginSpec().registry, behaviors: [] });
  await plugin.init();
  return plugin;
}

async function load(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: pdb, label: 'a' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  return (await plugin.builders.structure.createStructure(model)).data as Structure;
}

const atom = (s: Structure, i: number) =>
  StructureElement.Loci(s, [{ unit: s.units[0], indices: OrderedSet.ofSingleton(i as StructureElement.UnitIndex) }]);
const lines = { visualParams: { visuals: ['lines'] } } as const;

afterEach(() => vi.restoreAllMocks());

// The group is kept from cell events, not searched for per change (perf audit 2, 10/P5).
describe('measurement sync', () => {
  it('searches the tree for the measurement group at no change, a measurement added or removed', async () => {
    const plugin = await createPlugin();
    const s = await load(plugin);
    const m = plugin.managers.structure.measurement;

    const find = vi.spyOn(StateSelection, 'findTagInSubtree');
    await plugin.builders.data.rawData({ data: 'x', label: 'other' });
    expect(find).not.toHaveBeenCalled();
    expect(m.state.distances).toEqual([]);

    const d = await m.addDistance(atom(s, 0), atom(s, 1), lines);
    expect(m.state.distances.length).toBe(1);
    await load(plugin);
    expect(m.state.distances.length).toBe(1);
    await plugin.build().delete(d!.selection.ref).commit();
    expect(m.state.distances).toEqual([]);
    // the one search is `addDistance`'s own (`getGroup`, to put the measurement under the group)
    expect(find).toHaveBeenCalledTimes(1);
    plugin.dispose();
  });

  it('finds the group again once it is gone, and a group tagged in place', async () => {
    const plugin = await createPlugin();
    const s = await load(plugin);
    const m = plugin.managers.structure.measurement;

    await m.addDistance(atom(s, 0), atom(s, 1), lines);
    const group = StateSelection.findTagInSubtree(plugin.state.data.tree, StateTransform.RootRef, MeasurementGroupTag)!;
    await plugin.build().delete(group).commit();
    expect(m.state.distances).toEqual([]);

    // an untagged group, tagged by an update: its measurements are listed
    const other = plugin.build().toRoot().apply(CreateGroup, { label: 'Other' });
    await other.commit();
    await plugin.build().to(other.ref).tag(MeasurementGroupTag).commit();
    await m.addDistance(atom(s, 0), atom(s, 1), lines);
    expect(StateSelection.findTagInSubtree(plugin.state.data.tree, StateTransform.RootRef, MeasurementGroupTag)).toBe(
      other.ref,
    );
    expect(m.state.distances.length).toBe(1);
    plugin.dispose();
  });
});
