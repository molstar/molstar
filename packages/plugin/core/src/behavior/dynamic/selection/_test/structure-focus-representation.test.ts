/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateSelection } from '@molstar/core/state';
import { Structure } from '@molstar/model/model/structure';
import { PluginBehaviors } from '@molstar/plugin/behavior';
import {
  StructureFocusRepresentation,
  StructureFocusRepresentationTags,
} from '@molstar/plugin/behavior/dynamic/selection/structure-focus-representation';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';

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
ATOM 5 N N . GLY A 1 2 3.300 1.450 0.000 1 10 2 A 1
ATOM 6 C CA . GLY A 1 2 4.000 2.700 0.000 1 10 2 A 1
ATOM 7 C C . GLY A 1 2 5.500 2.600 0.000 1 10 2 A 1
ATOM 8 O O . GLY A 1 2 6.100 1.500 0.000 1 10 2 A 1
#
`;

async function createPlugin(behaviors: PluginSpec['behaviors']) {
  const plugin = new PluginContext({ registry: DefaultPluginSpec().registry, behaviors });
  await plugin.init();
  return plugin;
}

async function focusStructure(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: Cif, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  plugin.managers.structure.focus.setFromLoci(Structure.toStructureElementLoci(structure.data!));
  const select = (tag: StructureFocusRepresentationTags) =>
    plugin.state.data.select(StateSelection.Generators.root.subtree().withTag(tag));
  // the behavior builds the focus representations asynchronously
  for (let i = 0; i < 200 && select(StructureFocusRepresentationTags.SurrRepr).length === 0; i++) {
    await new Promise((resolve) => setTimeout(resolve, 25));
  }
  await new Promise((resolve) => setTimeout(resolve, 100));
  return select;
}

describe('structure focus representation behavior', () => {
  it('skips the interactions component when the Interactions representation is not registered', async () => {
    const plugin = await createPlugin([PluginSpec.Behavior(StructureFocusRepresentation)]);
    const select = await focusStructure(plugin);
    expect(select(StructureFocusRepresentationTags.TargetRepr)).toHaveLength(1);
    expect(select(StructureFocusRepresentationTags.SurrRepr)).toHaveLength(1);
    expect(select(StructureFocusRepresentationTags.SurrNciRepr)).toHaveLength(0);
    plugin.dispose();
  }, 30000);

  it('applies the interactions component when the Interactions behavior registers its representation', async () => {
    const plugin = await createPlugin([
      PluginSpec.Behavior(PluginBehaviors.CustomProps.Interactions),
      PluginSpec.Behavior(StructureFocusRepresentation),
    ]);
    const select = await focusStructure(plugin);
    expect(select(StructureFocusRepresentationTags.TargetRepr)).toHaveLength(1);
    expect(select(StructureFocusRepresentationTags.SurrNciRepr)).toHaveLength(1);
    plugin.dispose();
  }, 30000);
});
