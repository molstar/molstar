/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { CubeProvider } from '../../formats/volume/cube.js';
import { trackObjectRefs } from '../object-refs.js';

const cube = [
  'comment 1',
  'comment 2',
  '    1    0.000000    0.000000    0.000000',
  '    2    1.000000    0.000000    0.000000',
  '    2    0.000000    1.000000    0.000000',
  '    2    0.000000    0.000000    1.000000',
  '    1    0.000000    0.000000    0.000000    0.000000',
  '  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8',
  '',
].join('\n');

const cif = `data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 0.000 0.000 0.000 1 10 1 A 1
ATOM 2 C CA . GLY A 1 1 1.450 0.000 0.000 1 10 1 A 1
#
`;

async function createPlugin() {
  const plugin = new PluginContext({ behaviors: [], registry: DefaultPluginSpec().registry });
  await plugin.init();
  return plugin;
}

async function loadVolume(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: cube, label: 'blob.cube' });
  const parsed = await CubeProvider.parse(plugin, data);
  await CubeProvider.visuals!(plugin, parsed);
  return data;
}

async function loadStructure(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: cif, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  await plugin.builders.structure.createStructure(model);
  return data;
}

describe('object refs of a hierarchy manager', () => {
  it('follow objects created, updated, and removed', async () => {
    const plugin = await createPlugin();
    const refs = trackObjectRefs(plugin.managers.volume.hierarchy, plugin.state.data, (o) => SO.Volume.Data.is(o));
    await loadStructure(plugin);
    expect(refs.size).toBe(0);
    const v = await loadVolume(plugin);
    expect(refs.size).toBeGreaterThan(0);
    await plugin.managers.volume.hierarchy.remove([v.ref]);
    expect(refs.size).toBe(0);
    plugin.dispose();
  });

  it('keep the volume hierarchy through a volume added to and removed from a scene of structures', async () => {
    const plugin = await createPlugin();
    const volumes = () => plugin.managers.volume.hierarchy.current.volumes.map((v) => v.representations.length);
    const changes = vi.fn();
    plugin.managers.volume.hierarchy.behaviors.selection.subscribe(changes);
    changes.mockClear();

    await loadStructure(plugin);
    expect(volumes()).toEqual([]);
    expect(changes).not.toHaveBeenCalled();

    const v = await loadVolume(plugin);
    expect(volumes()).toEqual([1]);
    await loadStructure(plugin);
    expect(volumes()).toEqual([1]);

    await plugin.managers.volume.hierarchy.remove([v.ref]);
    expect(volumes()).toEqual([]);
    expect(plugin.managers.volume.hierarchy.current.refs.size).toBe(0);
    await loadStructure(plugin);
    expect(volumes()).toEqual([]);
    expect(plugin.managers.particles.hierarchy.current.lists).toEqual([]);
    plugin.dispose();
  });
});
