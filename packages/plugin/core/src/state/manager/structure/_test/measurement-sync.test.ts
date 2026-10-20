/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { StateSelection } from '@molstar/core/state';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';

const pdb = ['ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  '].join('\n');

async function createPlugin() {
  const plugin = new PluginContext({ registry: DefaultPluginSpec().registry, behaviors: [] });
  await plugin.init();
  return plugin;
}

describe('measurement sync', () => {
  it('searches the tree for the measurement group once a change', async () => {
    const plugin = await createPlugin();
    const data = await plugin.builders.data.rawData({ data: pdb, label: 'a' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    await plugin.builders.structure.createStructure(model);

    const find = vi.spyOn(StateSelection, 'findTagInSubtree');
    await plugin.builders.data.rawData({ data: 'x', label: 'other' });
    // One search for the change (six before, each a walk of the whole tree when there is no group).
    expect(find).toHaveBeenCalledTimes(1);
    expect(plugin.managers.structure.measurement.state.distances).toEqual([]);

    plugin.dispose();
  });
});
