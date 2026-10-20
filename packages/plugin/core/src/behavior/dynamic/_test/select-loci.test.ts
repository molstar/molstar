/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '@molstar/core/data/int';
import { MarkerAction } from '@molstar/core/util/marker-action';
import type { Representation } from '@molstar/graphics/repr/representation';
import { type Structure, StructureElement } from '@molstar/model/model/structure';
import { PluginBehaviors } from '@molstar/plugin/behavior';
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
#
`;

/** A canvas without WebGL: marks every representation added to it (canvas3d.ts `applyMark`), and records the marks */
function recordingCanvas() {
  const reprs = new Set<Representation.Any>();
  const marks: { action: MarkerAction; changed: boolean }[] = [];
  return {
    marks,
    canvas: {
      add: (r: Representation.Any) => reprs.add(r),
      remove: (r: Representation.Any) => reprs.delete(r),
      update() {},
      syncVisibility() {},
      requestDraw() {},
      commit() {},
      dispose() {},
      mark(reprLoci: Representation.Loci, action: MarkerAction) {
        let changed = false;
        for (const r of reprs) changed = r.mark(reprLoci.loci, action) || changed;
        marks.push({ action, changed });
      },
    },
  };
}

async function setup() {
  const plugin = new PluginContext({
    registry: DefaultPluginSpec().registry,
    behaviors: [PluginSpec.Behavior(PluginBehaviors.Representation.SelectLoci)],
  });
  await plugin.init();
  const c = recordingCanvas();
  (plugin as any).canvas3d = c.canvas;
  const data = await plugin.builders.data.rawData({ data: Cif, label: 'test' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { plugin, structure, marks: c.marks };
}

const atom = (s: Structure, i: number) =>
  StructureElement.Loci(s, [{ unit: s.units[0], indices: OrderedSet.ofSingleton(i as StructureElement.UnitIndex) }]);

describe('select loci behavior', () => {
  it('marks a new representation with the selection of its structure', async () => {
    const { plugin, structure, marks } = await setup();
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'ball-and-stick' });
    plugin.managers.interactivity.lociSelects.select({ loci: atom(structure.data!, 1) });
    marks.length = 0;
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
    // the new representation gets the selection (the one before has it already)
    expect(marks).toEqual([{ action: MarkerAction.Select, changed: true }]);
    plugin.dispose();
  });

  it('marks nothing for a new representation when its structure has nothing selected', async () => {
    const { plugin, structure, marks } = await setup();
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'ball-and-stick' });
    marks.length = 0;
    await plugin.builders.structure.representation.addRepresentation(structure, { type: 'spacefill' });
    // marking the empty selection on every representation changed nothing
    expect(marks).toEqual([]);
    plugin.dispose();
  });
});
