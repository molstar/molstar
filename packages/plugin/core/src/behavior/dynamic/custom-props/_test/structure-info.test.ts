/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Model, Structure } from '@molstar/model/model/structure';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { PluginSpec } from '@molstar/plugin/spec';
import { ModelFromTrajectory } from '../../../../state/transforms/structure/hierarchy.js';
import { StructureInfo } from '../structure-info.js';

const pdb = 'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  ';

async function createPlugin() {
  const plugin = new PluginContext({
    registry: DefaultPluginSpec().registry,
    behaviors: [PluginSpec.Behavior(StructureInfo)],
  });
  await plugin.init();
  return plugin;
}

async function load(plugin: PluginContext, label: string) {
  const data = await plugin.builders.data.rawData({ data: pdb, label });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { data, model: model.obj!.data, structure: structure.obj!.data };
}

describe('StructureInfo', () => {
  it('indexes models and structures and keeps their max index, lowered when a model goes', async () => {
    const plugin = await createPlugin();
    const loaded = [await load(plugin, 'a'), await load(plugin, 'b'), await load(plugin, 'c')];
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
    const d = await load(plugin, 'd');
    expect(Model.Index.get(d.model).value).toBe(2);
    expect(Model.MaxIndex.get(loaded[0].model).value).toBe(2);
    plugin.dispose();
  });

  it('gives a model and its structure replaced by another frame the max index at once, the structure its index kept', async () => {
    const plugin = await createPlugin();
    const twoModels = ['MODEL        1', pdb, 'ENDMDL', 'MODEL        2', pdb, 'ENDMDL'].join('\n');
    const data = await plugin.builders.data.rawData({ data: twoModels, label: 'nmr' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const modelCell = await plugin.builders.structure.createModel(trajectory);
    const structureCell = await plugin.builders.structure.createStructure(modelCell);
    await load(plugin, 'b');
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

  it('walks the tree only for models and structures, not for every object created', async () => {
    const plugin = await createPlugin();
    await load(plugin, 'a');
    const select = vi.spyOn(plugin.state.data, 'select');
    await plugin.builders.data.rawData({ data: 'x', label: 'other' });
    // The state's own select of the update, none of the behavior's (it made four before).
    expect(select).toHaveBeenCalledTimes(1);
    select.mockClear();
    await load(plugin, 'b');
    expect(select.mock.calls.length).toBeGreaterThan(4);
    plugin.dispose();
  });
});
