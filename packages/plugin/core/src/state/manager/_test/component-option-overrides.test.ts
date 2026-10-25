/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Material } from '@molstar/core/util/material';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { StructureComponentManager } from '../structure/component.js';

const Cif = `data_test
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
ATOM 3 C C . GLY A 1 1 2.000 1.400 0.000 1 10 1 A 1
ATOM 4 O O . GLY A 1 1 1.300 2.400 0.000 1 10 1 A 1
#
`;

const Glossy = { metalness: 0, roughness: 0.2, bumpiness: 0 };
const Shiny = { metalness: 0.5, roughness: 0.3, bumpiness: 0 };

// The component options are the representations' scene-wide defaults: a representation with a value of its own keeps
// it when they change (StructureComponentManager.setOptions).
async function setup() {
  const plugin = new PluginContext({ registry: DefaultPluginSpec().registry, behaviors: [] });
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: Cif, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  const component = await plugin.builders.structure.tryCreateComponentStatic(structure, 'all');
  const manager = plugin.managers.structure.component;
  const add = (type: string) =>
    plugin.builders.structure.representation.addRepresentation(component!, {
      type: plugin.representation.structure.registry.get(type),
      typeParams: { material: manager.state.options.materialStyle, quality: manager.state.options.visualQuality },
    });
  const a = await add('ball-and-stick');
  const b = await add('spacefill');
  const params = (ref: string) => plugin.state.data.cells.get(ref)!.transform.params!.type.params;
  return { plugin, manager, a: a.ref, b: b.ref, params };
}

describe('StructureComponentManager option overrides', () => {
  it('moves the representations that follow a changed option and keeps those with their own value', async () => {
    const { plugin, manager, a, b, params } = await setup();
    // b gets its own material
    await plugin
      .build()
      .to(b)
      .update(StructureRepresentation3D, (old) => {
        old.type.params.material = Shiny;
      })
      .commit();

    await manager.setOptions({ ...manager.state.options, materialStyle: Glossy });
    expect(Material.areEqual(params(a).material, Glossy)).toBe(true);
    expect(Material.areEqual(params(b).material, Shiny)).toBe(true);

    // another option changes: b's material stays its own, both take the quality
    await manager.setOptions({ ...manager.state.options, visualQuality: 'high' });
    expect(params(a).quality).toBe('high');
    expect(params(b).quality).toBe('high');
    expect(Material.areEqual(params(b).material, Shiny)).toBe(true);
    plugin.dispose();
  });

  it('writes the options in reset to every representation, overrides included', async () => {
    const { plugin, manager, b, params } = await setup();
    await plugin
      .build()
      .to(b)
      .update(StructureRepresentation3D, (old) => {
        old.type.params.material = Shiny;
      })
      .commit();

    await manager.setOptions(manager.state.options, { reset: ['materialStyle'] });
    expect(Material.areEqual(params(b).material, manager.state.options.materialStyle)).toBe(true);
    plugin.dispose();
  });

  it('reads overrides and values in the options terms', async () => {
    const { plugin, manager, a, b, params } = await setup();
    await plugin
      .build()
      .to(b)
      .update(StructureRepresentation3D, (old) => {
        old.type.params.material = Shiny;
        old.type.params.ignoreHydrogens = true;
        old.type.params.ignoreHydrogensVariant = 'non-polar';
      })
      .commit();
    const options = manager.state.options;
    expect(StructureComponentManager.reprOptionOverrides('ball-and-stick', params(a), options)).toEqual([]);
    expect(StructureComponentManager.reprOptionOverrides('spacefill', params(b), options)).toEqual([
      'hydrogens',
      'materialStyle',
    ]);
    expect(StructureComponentManager.reprOptionValue('spacefill', params(b), 'hydrogens')).toBe('only-polar');
    // a cartoon takes no solid surface
    expect(StructureComponentManager.reprHasOption('cartoon', params(a), 'solidSurface', options)).toBeUndefined();
    plugin.dispose();
  });
});
