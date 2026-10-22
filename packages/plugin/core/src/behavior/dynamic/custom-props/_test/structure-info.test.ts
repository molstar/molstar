/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { StateSelection } from '@molstar/core/state';
import { Model, Structure } from '@molstar/model/model/structure';
import { PluginBehaviors } from '@molstar/plugin/behavior';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import { StructureInfo } from '../structure-info.js';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { RawData } from '@molstar/plugin/state/transforms/data/fetch';
import { ModelFromTrajectory, StructureFromModel } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { ParseCif } from '../../../../state/formats/cif.js';
import { TrajectoryFromMmCif } from '../../../../state/formats/trajectory/mmcif.js';

const pdb = 'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  ';

async function createPluginWithBehavior() {
  const plugin = new PluginContext({
    registry: DefaultPluginSpec().registry,
    behaviors: [PluginSpec.Behavior(StructureInfo)],
  });
  await plugin.init();
  return plugin;
}

async function loadPdb(plugin: PluginContext, label: string) {
  const data = await plugin.builders.data.rawData({ data: pdb, label });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { data, model: model.obj!.data, structure: structure.obj!.data };
}

// the spies are on shared module objects
afterEach(() => vi.restoreAllMocks());

describe('StructureInfo', () => {
  it('indexes models and structures and keeps their max index, lowered when a model goes', async () => {
    const plugin = await createPluginWithBehavior();
    const loaded = [await loadPdb(plugin, 'a'), await loadPdb(plugin, 'b'), await loadPdb(plugin, 'c')];
    expect(loaded.map((l) => Model.Index.get(l.model).value)).toEqual([0, 1, 2]);
    expect(loaded.map((l) => Structure.Index.get(l.structure).value)).toEqual([0, 1, 2]);
    for (const l of loaded) {
      expect(Model.MaxIndex.get(l.model).value).toBe(2);
      expect(Structure.MaxIndex.get(l.structure).value).toBe(2);
    }
    await plugin.build().delete(loaded[2].data).commit();
    for (const l of loaded.slice(0, 2)) {
      expect(Model.MaxIndex.get(l.model).value).toBe(1);
      expect(Structure.MaxIndex.get(l.structure).value).toBe(1);
    }
    // A new model takes the next index after the highest left.
    const d = await loadPdb(plugin, 'd');
    expect(Model.Index.get(d.model).value).toBe(2);
    expect(Model.MaxIndex.get(loaded[0].model).value).toBe(2);
    plugin.dispose();
  });

  it('gives a model and its structure replaced by another frame the max index at once, the structure its index kept', async () => {
    const plugin = await createPluginWithBehavior();
    const twoModels = ['MODEL        1', pdb, 'ENDMDL', 'MODEL        2', pdb, 'ENDMDL'].join('\n');
    const data = await plugin.builders.data.rawData({ data: twoModels, label: 'nmr' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const modelCell = await plugin.builders.structure.createModel(trajectory);
    const structureCell = await plugin.builders.structure.createStructure(modelCell);
    await loadPdb(plugin, 'b');
    const first = modelCell.obj!.data;
    const firstStructure = structureCell.obj!.data;
    expect(Structure.Index.get(firstStructure).value).toBe(0);
    expect(Model.MaxIndex.get(first).value).toBe(1);
    await plugin
      .build()
      .to(modelCell.ref)
      .update(ModelFromTrajectory, (old: any) => ({ ...old, modelIndex: 1 }))
      .commit();
    const model = plugin.state.data.cells.get(modelCell.ref)!.obj!.data as Model;
    expect(model).not.toBe(first);
    expect(Model.MaxIndex.get(model).value).toBe(1);
    // Its structure, remapped onto the new model in place, keeps its index and has the max too.
    const structure = plugin.state.data.cells.get(structureCell.ref)!.obj!.data as Structure;
    expect(structure).not.toBe(firstStructure);
    expect(Structure.Index.get(structure).value).toBe(0);
    expect(Structure.MaxIndex.get(structure).value).toBe(1);
    plugin.dispose();
  });

  it('walks the tree for no object created, models and structures included', async () => {
    const plugin = await createPluginWithBehavior();
    await loadPdb(plugin, 'a');
    const select = vi.spyOn(plugin.state.data, 'select');
    await plugin.builders.data.rawData({ data: 'x', label: 'other' });
    // The state's own select of the update, none of the behavior's (it made four before).
    expect(select).toHaveBeenCalledTimes(1);
    // nor a walk for a model or a structure: the behavior keeps its roots (14/F6)
    const roots = vi.spyOn(StateSelection.Generators, 'rootsOfType');
    await loadPdb(plugin, 'b');
    expect(roots).not.toHaveBeenCalled();
    plugin.dispose();
  });
});

/** An mmCIF entry of `chains` chains of two glycines each */
function cif(name: string, chains: number) {
  const rows: string[] = [];
  let id = 1;
  for (let c = 0; c < chains; c++) {
    const asym = String.fromCharCode(65 + c);
    for (let r = 1; r <= 2; r++) {
      for (const [atom, el] of [
        ['N', 'N'],
        ['CA', 'C'],
        ['C', 'C'],
        ['O', 'O'],
      ]) {
        rows.push(`ATOM ${id} ${el} ${atom} . GLY ${asym} 1 ${r} ${id * 1.1} ${c * 5} 0 1 10 ${r} ${asym} 1`);
        id++;
      }
    }
  }
  return `data_${name}
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
${rows.join('\n')}
#
`;
}

async function createPlugin() {
  const plugin = new PluginContext({
    registry: DefaultPluginSpec().registry,
    behaviors: [PluginSpec.Behavior(PluginBehaviors.CustomProps.StructureInfo)],
  });
  await plugin.init();
  return plugin;
}

async function load(plugin: PluginContext, name: string, chains: number) {
  const data = await plugin.builders.data.rawData({ data: cif(name, chains), label: name });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  const model = await plugin.builders.structure.createModel(trajectory);
  const props = await plugin.builders.structure.insertModelProperties(model);
  const structure = await plugin.builders.structure.createStructure(props);
  return { data, structure };
}

/** The structure info of every model and structure object in the state, in cell order */
function info(plugin: PluginContext) {
  const rows: string[] = [];
  for (const cell of plugin.state.data.cells.values()) {
    const label = cell.obj?.label;
    if (SO.Molecule.Model.is(cell.obj)) {
      const m = cell.obj.data;
      const o = Model.AsymIdOffset.get(m).value;
      rows.push(
        `model ${label}: index=${Model.Index.get(m).value} max=${Model.MaxIndex.get(m).value} asym=${o ? `${o.auth}/${o.label}` : o}`,
      );
    } else if (SO.Molecule.Structure.is(cell.obj)) {
      const s = cell.obj.data;
      rows.push(`structure ${label}: index=${Structure.Index.get(s).value} max=${Structure.MaxIndex.get(s).value}`);
    }
  }
  return rows;
}

describe('structure info behavior', () => {
  it('numbers models and structures and keeps their max index over loads, updates and removals', async () => {
    const plugin = await createPlugin();
    const steps: Record<string, string[]> = {};

    const a = await load(plugin, 'a', 2);
    const b = await load(plugin, 'b', 3);
    const c = await load(plugin, 'c', 1);
    steps.loaded = info(plugin);

    // a representation is created: an object that is not a model or a structure
    await plugin.builders.structure.representation.addRepresentation(a.structure, { type: 'ball-and-stick' });
    steps.representation = info(plugin);

    // the structure is updated (recreated): its index carries over
    await plugin
      .build()
      .to(b.structure)
      .update({ type: { name: 'model', params: {} } })
      .commit();
    steps.updated = info(plugin);

    // the last model goes away; the max index is lowered at once (#2008)
    await plugin.managers.structure.hierarchy.remove([c.data.ref]);
    steps.removed = info(plugin);
    await plugin.builders.structure.representation.addRepresentation(b.structure, { type: 'cartoon' });
    steps.createdAfterRemove = info(plugin);

    await load(plugin, 'd', 2);
    steps.loadedAgain = info(plugin);

    expect(steps).toMatchSnapshot();
    plugin.dispose();
  });

  it('numbers models and structures created in one update as those created one by one', async () => {
    const plugin = await createPlugin();
    await load(plugin, 'a', 1);
    const b = plugin.build();
    for (const [name, chains] of [
      ['b', 2],
      ['c', 3],
    ] as const) {
      b.toRoot()
        .apply(RawData, { data: cif(name, chains), label: name })
        .apply(ParseCif)
        .apply(TrajectoryFromMmCif)
        .apply(ModelFromTrajectory, { modelIndex: 0 })
        .apply(StructureFromModel);
    }
    await b.commit();
    expect(info(plugin)).toMatchSnapshot();
    plugin.dispose();
  });

  it('takes the models and structures loaded before it was registered', async () => {
    const plugin = new PluginContext({ registry: DefaultPluginSpec().registry, behaviors: [] });
    await plugin.init();
    await load(plugin, 'a', 2);
    await load(plugin, 'b', 1);
    await plugin.state.updateBehavior(PluginBehaviors.CustomProps.StructureInfo, () => {});
    await load(plugin, 'c', 3);
    // a component: not a root structure
    const s = plugin.managers.structure.hierarchy.current.structures[0].cell;
    await plugin.builders.structure.tryCreateComponentStatic(s, 'polymer');
    expect(info(plugin)).toMatchSnapshot();
    plugin.dispose();
  });

  it('walks the state for no created object once registered, components and representations included', async () => {
    const plugin = await createPlugin();
    const a = await load(plugin, 'a', 2);
    await load(plugin, 'b', 1);
    const select = vi.spyOn(StateSelection.Generators, 'rootsOfType');
    await plugin.builders.structure.representation.addRepresentation(a.structure, { type: 'ball-and-stick' });
    // a component is a structure too, but not a root one
    const polymer = await plugin.builders.structure.tryCreateComponentStatic(a.structure, 'polymer');
    await plugin.builders.structure.representation.addRepresentation(polymer!, { type: 'cartoon' });
    expect(select).not.toHaveBeenCalled();
    plugin.dispose();
  });

  it('visits the roots for a root model only: a component or a root structure visits none (perf audit 2, 10/P4)', async () => {
    const plugin = await createPlugin();
    const a = await load(plugin, 'a', 2);
    await load(plugin, 'b', 1);
    await load(plugin, 'c', 3);
    const index = vi.spyOn(Model.Index, 'get');
    const modelMax = vi.spyOn(Model.MaxIndex, 'get');
    const structureMax = vi.spyOn(Structure.MaxIndex, 'get');
    const polymer = await plugin.builders.structure.tryCreateComponentStatic(a.structure, 'polymer');
    await plugin.builders.structure.representation.addRepresentation(polymer!, { type: 'cartoon' });
    await plugin.build().delete(polymer!).commit();
    expect([index, modelMax, structureMax].map((s) => s.mock.calls.length)).toEqual([0, 0, 0]);
    // a second root structure of a model: its own max, no other root's
    await plugin.builders.structure.createStructure(a.structure.cell!.transform.parent);
    expect([index, modelMax, structureMax].map((s) => s.mock.calls.length)).toEqual([0, 0, 1]);
    plugin.dispose();
  });
});
