/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { StructureMeasurementManagerState } from '../measurement.js';

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

async function createPlugin() {
  const plugin = new PluginContext({ registry: DefaultPluginSpec().registry, behaviors: [] });
  await plugin.init();
  return plugin;
}

async function loadStructure(plugin: PluginContext, label: string) {
  const data = await plugin.builders.data.rawData({ data: Cif, label });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  return (await plugin.builders.structure.createStructure(model)).data as Structure;
}

const atom = (s: Structure, i: number) =>
  StructureElement.Loci(s, [{ unit: s.units[0], indices: OrderedSet.ofSingleton(i as StructureElement.UnitIndex) }]);

/** the measurement state as labels of its cells (refs are random) */
function summary(state: StructureMeasurementManagerState) {
  const labels = (cells: readonly any[]) => cells.map((c) => `${c.transform.transformer.id}:${c.obj?.label}`);
  return {
    labels: labels(state.labels),
    distances: labels(state.distances),
    angles: labels(state.angles),
    dihedrals: labels(state.dihedrals),
    orientations: labels(state.orientations),
    planes: labels(state.planes),
  };
}

describe('structure measurement manager', () => {
  it('lists the measurements of the measurement group after every change', async () => {
    const plugin = await createPlugin();
    const s = await loadStructure(plugin, 'a');
    const t = await loadStructure(plugin, 'b');
    const m = plugin.managers.structure.measurement;
    const steps: Record<string, unknown> = { empty: summary(m.state) };

    // no text visuals: their font atlas needs a Canvas module in Node (a label fails, so it is listed nowhere)
    const d1 = await m.addDistance(atom(s, 0), atom(s, 1), { visualParams: { visuals: ['lines'] } });
    steps.distance = summary(m.state);
    await m.addDistance(atom(s, 2), atom(t, 3), { visualParams: { visuals: ['lines'] } });
    await m.addLabel(atom(t, 0));
    await m.addAngle(atom(s, 0), atom(s, 1), atom(s, 2), { visualParams: { visuals: ['vectors', 'sector'] } });
    await m.addDihedral(atom(s, 0), atom(s, 1), atom(s, 2), atom(s, 3), {
      visualParams: { visuals: ['extenders', 'arms', 'sector'] },
    });
    await m.addOrientation([atom(s, 0)]);
    await m.addPlane([atom(t, 1)]);
    steps.all = summary(m.state);

    // an unrelated change: another structure
    await loadStructure(plugin, 'c');
    steps.unrelated = summary(m.state);

    await plugin.build().delete(d1!.selection.ref).commit();
    steps.removed = summary(m.state);

    expect(steps).toMatchSnapshot();
    plugin.dispose();
  });
});
