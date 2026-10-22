/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { type Model, Structure, StructureElement } from '@molstar/model/model/structure';
import { OrderedSet } from '@molstar/core/data/int';
import { CustomModelProperty } from '@molstar/model/props/common/custom-model-property';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginContext } from '@molstar/plugin/context';
import { ModelFromTrajectory } from '@molstar/plugin/state/transforms/structure/hierarchy';
import { Color } from '@molstar/core/util/color';
import { setStructureOverpaint } from '@molstar/plugin/state/helpers/structure-overpaint';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import { createPlaybackPlugin, pdbEnsemble, play } from './harness.js';

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

describe('the components of a trajectory frame', () => {
  async function components() {
    const env = await createPlaybackPlugin(DefaultPluginSpec());
    const { plugin } = env;
    const { model } = await loadEnsemble(plugin);
    const structure = await plugin.builders.structure.createStructure(model);
    const polymer = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
    const water = await plugin.builders.structure.tryCreateComponentStatic(structure, 'water');
    const data = (ref: string) => plugin.state.data.cells.get(ref)!.obj!.data as Structure;
    return {
      ...env,
      model,
      root: () => data(structure.ref),
      polymer: () => data(polymer!.ref),
      water: () => data(water!.ref),
    };
  }

  it('share the remapped root (T2)', async () => {
    const { plugin, model, root, polymer, water } = await components();
    expect(polymer().root).toBe(root());
    for (const i of [1, 2]) {
      await setFrame(plugin, model.ref, i);
      expect(root().model.modelNum).toBe(i + 1);
      expect(polymer().root).toBe(root());
      expect(water().root).toBe(root());
      expect(Structure.areRootsEquivalent(polymer(), water())).toBe(true);
      // the root's property data, which components inherit (secondary structure, ASA)
      expect(polymer().inheritedPropertyData).toBe(root().inheritedPropertyData);
    }
  });

  it("take the frame's coordinates", async () => {
    const { plugin, model, root, polymer } = await components();
    await setFrame(plugin, model.ref, 2);
    const u = polymer().units[0];
    const r = root().unitMap.get(u.id);
    expect(u.conformation.x(u.elements[0])).toBe(r.conformation.x(u.elements[0]));
    expect(u.model).toBe(root().model);
  });
});

describe('a trajectory frame with the default preset and an overpaint', () => {
  it('updates every cell in place (T3)', async () => {
    const { plugin, counts } = await createPlaybackPlugin(DefaultPluginSpec());
    const data = await plugin.builders.data.rawData({ data: pdbEnsemble(3, 6, true), label: 'ensemble' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    const hierarchy = () => plugin.managers.structure.hierarchy.current;
    await setStructureOverpaint(plugin, hierarchy().structures[0].components, Color(0xff0000), async (s) =>
      StructureElement.Loci(s, [{ unit: s.units[0], indices: OrderedSet.ofBounds(0, 5) as any }]),
    );
    const recreated = counts.recreated;
    const created = counts.created;
    await setFrame(plugin, hierarchy().models[0].cell.transform.ref, 2);
    expect(counts.recreated - recreated).toBe(0);
    expect(counts.created - created).toBe(0);
    const root = hierarchy().structures[0].cell.obj!.data;
    expect(root.model.modelNum).toBe(3);
    for (const c of hierarchy().structures[0].components) {
      expect(c.cell.obj!.data.root).toBe(root);
      for (const o of c.representations[0].overpaint ? [c.representations[0].overpaint] : []) {
        const overpaint = (o.cell.obj!.data.state as any).overpaint;
        expect(overpaint.layers[0].loci.structure.root).toBe(root);
      }
    }
  });
});

describe('a focused trajectory frame', () => {
  async function focused() {
    const env = await createPlaybackPlugin(DefaultPluginSpec());
    const { plugin } = env;
    const data = await plugin.builders.data.rawData({ data: pdbEnsemble(3, 6, true), label: 'ensemble' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    const s = plugin.managers.structure.hierarchy.current.structures[0].cell.obj!.data;
    const unit = s.units[0];
    plugin.managers.structure.focus.setFromLoci(
      StructureElement.Loci(s, [{ unit, indices: OrderedSet.ofBounds(10, 15) as any }]),
    );
    await flush(plugin);
    const model = plugin.managers.structure.hierarchy.current.models[0].cell.transform.ref;
    return { ...env, model };
  }

  it('commits one update, not a second one for the focus (T5)', async () => {
    const { plugin, model } = await focused();
    const updates = vi.fn();
    plugin.state.data.events.changed.subscribe(updates);
    await setFrame(plugin, model, 1);
    await flush(plugin);
    expect(updates).toHaveBeenCalledTimes(1);
    // the target selection followed the frame in that update
    const target = [...plugin.state.data.cells.values()].find((c) => c.obj?.label === '[Focus] Target')!;
    expect((target.obj!.data as Structure).model.modelNum).toBe(2);
  });

  it('can be disposed right after a frame (B1)', async () => {
    const { plugin, model } = await focused();
    const errors = vi.spyOn(console, 'error').mockImplementation(() => {});
    await setFrame(plugin, model, 1);
    plugin.dispose();
    await new Promise((r) => setTimeout(r, 20));
    expect(errors.mock.calls.map((c) => String(c[0]))).toEqual([]);
  });
});

describe('a playing trajectory', () => {
  it('commits an update only on the ticks whose frame changes (T5)', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    await loadEnsemble(plugin, 3);
    const updates = vi.fn();
    plugin.state.data.events.changed.subscribe(updates);
    const busy = vi.fn();
    plugin.behaviors.state.isUpdating.subscribe(busy);
    busy.mockClear();
    // 3 frames over 1 s at 60 Hz: 60 ticks (the last at 983 ms), the frame changes twice (0 -> 1 -> 2)
    await plugin.managers.animation.play(AnimateModelIndex, {
      mode: { name: 'loop', params: { direction: 'forward' } },
      duration: { name: 'fixed', params: { durationInS: 1 } },
    } as any);
    const warn = vi.spyOn(console, 'warn');
    await play(plugin, 1000);
    expect(updates).toHaveBeenCalledTimes(2);
    expect(busy.mock.calls.filter((c) => c[0] === true)).toHaveLength(2);
    expect(warn).not.toHaveBeenCalled();
  });
});

/** until the state has no update queued */
async function flush(plugin: PluginContext) {
  for (let i = 0; i < 10; i++) {
    await new Promise((r) => setTimeout(r, 0));
    if (!plugin.state.data.behaviors.isUpdating.value) return;
  }
}
