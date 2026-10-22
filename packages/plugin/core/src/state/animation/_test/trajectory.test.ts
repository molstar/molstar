/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import type { Model } from '@molstar/model/model/structure';
import { CustomModelProperty } from '@molstar/model/props/common/custom-model-property';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginContext } from '@molstar/plugin/context';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { createPlaybackPlugin, pdbEnsemble } from './harness.js';

async function loadEnsemble(plugin: PluginContext, models = 3, residues = 4) {
  const data = await plugin.builders.data.rawData({ data: pdbEnsemble(models, residues, true), label: 'ensemble' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  return { trajectory, model };
}

const setFrame = (plugin: PluginContext, ref: string, modelIndex: number) =>
  plugin
    .build()
    .to(ref)
    .update(ModelFromTrajectory, () => ({ modelIndex }))
    .commit();

/** A static model property whose value holds an asset */
function propertyWithAsset() {
  const dispose = vi.fn();
  const provider = CustomModelProperty.createProvider({
    label: 'With asset',
    descriptor: CustomPropertyDescriptor({ name: `with-asset-${Math.random()}` }),
    type: 'static',
    defaultParams: {},
    getParams: () => ({}),
    isApplicable: () => true,
    obtain: async () => ({ value: 1, assets: [{ dispose }] }),
  });
  return { provider, dispose };
}

describe('a trajectory frame', () => {
  it('keeps the custom properties the frames share (T1)', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const { provider, dispose } = propertyWithAsset();
    const first = model.obj!.data as Model;
    await provider.attach({ runtime: undefined as any, assetManager: plugin.managers.asset }, first);
    for (const i of [1, 2, 0]) await setFrame(plugin, model.ref, i);
    const now = plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    expect(now.customProperties).toBe(first.customProperties);
    expect(provider.get(now).value).toBe(1);
    expect(dispose).not.toHaveBeenCalled();
    // deleting the model cell still disposes them
    await plugin.build().delete(model.ref).commit();
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('updates the model cell in place', async () => {
    const { plugin, counts } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const cell = plugin.state.data.cells.get(model.ref)!;
    const obj = cell.obj;
    const recreated = counts.recreated;
    await setFrame(plugin, model.ref, 2);
    expect((cell.obj!.data as Model).modelNum).toBe(3);
    expect(cell.obj!.label).toBe('Model 3');
    expect(cell.obj).toBe(obj);
    expect(counts.recreated - recreated).toBe(0);
  });
});
