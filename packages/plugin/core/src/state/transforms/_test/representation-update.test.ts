/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What an update of a structure representation does again, as an animation frame updates it: a frame of an MVS
// story changes a representation's own params (opacity, emissive, clipping), a trajectory frame its structure.

import { describe, expect, it, vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { pdbEnsemble } from '../../animation/_test/harness.js';

async function representation() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: pdbEnsemble(2, 6), label: 'ensemble' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  const repr = await plugin.builders.structure.representation.addRepresentation(structure, {
    type: 'ball-and-stick',
    color: 'chain-id',
  });
  const colorThemes = plugin.representation.structure.themes.colorThemeRegistry;
  const create = vi.spyOn(colorThemes, 'create');
  const update = (f: (old: any) => any) => plugin.build().to(repr).update(StructureRepresentation3D, f).commit();
  return { plugin, model, repr, create, update };
}

describe('a structure representation update', () => {
  it('keeps its themes when only its own params change (R1)', async () => {
    const { plugin, repr, create, update } = await representation();
    await update((old) => ({ ...old, type: { ...old.type, params: { ...old.type.params, alpha: 0.5 } } }));
    expect(create).not.toHaveBeenCalled();
    const obj = plugin.state.data.cells.get(repr.ref)!.obj!;
    expect(obj.data.repr.props.alpha).toBe(0.5);
  });

  it('makes the themes again when their params change', async () => {
    const { create, update } = await representation();
    await update((old) => ({ ...old, colorTheme: { name: 'element-symbol', params: {} } }));
    expect(create).toHaveBeenCalledTimes(1);
  });

  it('makes the themes again on a new structure', async () => {
    const { plugin, model, create } = await representation();
    await plugin
      .build()
      .to(model)
      .update(ModelFromTrajectory, () => ({ modelIndex: 1 }))
      .commit();
    expect(create).toHaveBeenCalledTimes(1);
  });
});
