/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { CustomPropertyDescriptor } from '@molstar/model/model/custom-property';
import { type Model, Structure, StructureElement, Unit } from '@molstar/model/model/structure';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { OrderedSet } from '@molstar/core/data/int';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { CustomModelProperty } from '@molstar/model/props/common/custom-model-property';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginContext } from '@molstar/plugin/context';
import {
  ModelFromTrajectory,
  ModelWithCoordinates,
  TransformStructureConformation,
} from '@molstar/plugin/state/transforms/structure/hierarchy';
import { Time } from '@molstar/model/model/structure/coordinates';
import { Color } from '@molstar/core/util/color';
import { setStructureOverpaint } from '@molstar/plugin/state/helpers/structure-overpaint';
import { OverpaintStructureRepresentation3DFromBundle } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import type { Overpaint } from '@molstar/graphics/theme/overpaint';
import { AnimateModelIndex } from '@molstar/plugin/state/animation/built-in/model-index';
import { RawData } from '@molstar/plugin/state/transforms/data/fetch';
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

/** A model property whose value holds an asset: static (the frames share it) or dynamic (each frame its own) */
function propertyWithAsset(type: 'static' | 'dynamic' = 'static') {
  const dispose = vi.fn();
  const provider = CustomModelProperty.createProvider({
    label: 'With asset',
    descriptor: CustomPropertyDescriptor({ name: `with-asset-${Math.random()}` }),
    type,
    defaultParams: {},
    getParams: () => ({}),
    isApplicable: () => true,
    obtain: async () => ({ value: 1, assets: [{ dispose }] }),
  });
  return { provider, dispose };
}

describe('a trajectory frame', () => {
  it('keeps the custom properties the frames share', async () => {
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

  it("disposes the old frame's own (dynamic) property assets", async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const { provider, dispose } = propertyWithAsset('dynamic');
    await provider.attach({ runtime: undefined as any, assetManager: plugin.managers.asset }, model.obj!.data as Model);
    await setFrame(plugin, model.ref, 1);
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('disposes the properties of a frame that does not share them (another topology)', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    // two models of 4 and 5 residues: each has its own custom properties
    const text = [pdbEnsemble(1, 4), pdbEnsemble(1, 5).replace('MODEL        1', 'MODEL        2')]
      .map((t) => t.replace(/\nEND$/, ''))
      .join('\n');
    const data = await plugin.builders.data.rawData({ data: `${text}\nEND`, label: 'two topologies' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const first = model.obj!.data as Model;
    const { provider, dispose } = propertyWithAsset();
    await provider.attach({ runtime: undefined as any, assetManager: plugin.managers.asset }, first);
    await setFrame(plugin, model.ref, 1);
    const second = plugin.state.data.cells.get(model.ref)!.obj!.data as Model;
    expect(second.customProperties).not.toBe(first.customProperties);
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('disposes them when the trajectory is a new one (the data read again)', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const data = await plugin.builders.data.rawData({ data: pdbEnsemble(3, 4, true), label: 'ensemble' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const { provider, dispose } = propertyWithAsset();
    await provider.attach({ runtime: undefined as any, assetManager: plugin.managers.asset }, model.obj!.data as Model);
    await plugin
      .build()
      .to(data)
      .update(RawData, (old) => ({ ...old, data: pdbEnsemble(3, 5, true) }))
      .commit();
    expect((plugin.state.data.cells.get(model.ref)!.obj!.data as Model).atomicHierarchy.residues._rowCount).toBe(6);
    expect(dispose).toHaveBeenCalledTimes(1);
  });

  it('reaches a model-with-coordinates decorator below it', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const m = model.obj!.data as Model;
    const n = m.atomicHierarchy.atoms._rowCount;
    const frame = {
      elementCount: n,
      time: Time(0, 'step'),
      x: new Float32Array(n).fill(1),
      y: new Float32Array(n),
      z: new Float32Array(n),
      xyzOrdering: { isIdentity: true },
    };
    const withCoordinates = await plugin
      .build()
      .to(model)
      .apply(ModelWithCoordinates, { atomicCoordinateFrame: frame, frameIndex: 0, frameCount: 1 })
      .commit();
    await setFrame(plugin, model.ref, 2);
    const decorated = plugin.state.data.cells.get(withCoordinates.ref)!.obj!.data as Model;
    expect(decorated.modelNum).toBe(3);
    expect(decorated.atomicConformation.x[0]).toBe(1);
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

  it('share the remapped root', async () => {
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

  it('run the query again with bonds computed per frame (dynamic bonds)', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const structure = await plugin.builders.structure.createStructure(model, {
      name: 'model',
      params: { dynamicBonds: true },
    } as any);
    // any static component (the "plus connected" queries, ligand and branched, follow the bonds)
    const water = await plugin.builders.structure.tryCreateComponentStatic(structure, 'water');
    const recreated: string[] = [];
    plugin.state.data.events.object.updated.subscribe((e) => e.action === 'recreate' && recreated.push(e.ref));
    await setFrame(plugin, model.ref, 1);
    expect(recreated).toContain(water!.ref);
  });

  it("share the root's units where they are whole", async () => {
    const { plugin, model, root, polymer } = await components();
    expect(polymer().units[0]).toBe(root().unitMap.get(polymer().units[0].id));
    await setFrame(plugin, model.ref, 1);
    const u = polymer().units[0];
    // a whole unit of the root is the root's unit: its boundary and lookup are computed once per frame
    expect(u).toBe(root().unitMap.get(u.id));
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

describe('a trajectory frame that also moves the structure', () => {
  it('puts the components where the structure is', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const structure = await plugin.builders.structure.createStructure(model);
    const moved = await plugin
      .build()
      .to(structure)
      .apply(TransformStructureConformation, {
        transform: { name: 'matrix', params: { data: Mat4.identity(), transpose: false } },
      })
      .commit();
    const polymer = await plugin.builders.structure.tryCreateComponentStatic(moved, 'polymer');
    // one update: the next frame, and a translation of the structure
    const b = plugin.build();
    b.to(model).update(ModelFromTrajectory, () => ({ modelIndex: 1 }));
    b.to(moved).update(TransformStructureConformation, () => ({
      transform: {
        name: 'matrix',
        params: { data: Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0)), transpose: false },
      },
    }));
    await b.commit();
    const root = plugin.state.data.cells.get(moved.ref)!.obj!.data as Structure;
    const component = plugin.state.data.cells.get(polymer!.ref)!.obj!.data as Structure;
    const u = component.units[0];
    expect(Unit.areOperatorsEqual(u, root.unitMap.get(u.id))).toBe(true);
    expect(u.conformation.x(u.elements[0])).toBeCloseTo(root.unitMap.get(u.id).conformation.x(u.elements[0]), 5);
    expect(u.model).toBe(root.model);
  });

  it('keeps a whole-unit component of a moved structure where the structure is (perf-audit-2 05/B9)', async () => {
    const { plugin } = await createPlaybackPlugin(DefaultPluginSpec());
    const { model } = await loadEnsemble(plugin);
    const structure = await plugin.builders.structure.createStructure(model);
    const moved = await plugin
      .build()
      .to(structure)
      .apply(TransformStructureConformation, {
        transform: {
          name: 'matrix',
          params: { data: Mat4.fromTranslation(Mat4(), Vec3.create(10, 0, 0)), transpose: false },
        },
      })
      .commit();
    const polymer = await plugin.builders.structure.tryCreateComponentStatic(moved, 'polymer');
    // only the frame changes: the component is remapped, not recreated
    await setFrame(plugin, model.ref, 1);
    const root = plugin.state.data.cells.get(moved.ref)!.obj!.data as Structure;
    const component = plugin.state.data.cells.get(polymer!.ref)!.obj!.data as Structure;
    expect(component.model.modelNum).toBe(root.model.modelNum);
    const u = component.units[0];
    const r = root.unitMap.get(u.id);
    expect(Unit.areOperatorsEqual(u, r)).toBe(true);
    expect(u.conformation.x(u.elements[0])).toBeCloseTo(r.conformation.x(u.elements[0]), 5);
  });
});

describe('a trajectory frame with the default preset and an overpaint', () => {
  it('updates every cell in place', async () => {
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
    for (const c of hierarchy().structures[0].components) expect(c.cell.obj!.data.root).toBe(root);
    // the overpaint layers are on the new frame's structure (updated in place on it: the effects' structure change)
    const overpaints = [...plugin.state.data.cells.values()].filter(
      (c) => c.transform.transformer === OverpaintStructureRepresentation3DFromBundle,
    );
    expect(overpaints.length).toBeGreaterThan(0);
    for (const o of overpaints) {
      const overpaint = (o.obj!.data.state as { overpaint: Overpaint }).overpaint;
      expect(overpaint.layers).toHaveLength(1);
      expect((overpaint.layers[0].loci as StructureElement.Loci).structure.root).toBe(root);
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

  it('commits one update, not a second one for the focus', async () => {
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

  it('can be disposed right after a frame', async () => {
    const { plugin, model } = await focused();
    const errors = vi.spyOn(console, 'error').mockImplementation(() => {});
    await setFrame(plugin, model, 1);
    plugin.dispose();
    // the update the focus queued runs after the dispose, if any
    for (let i = 0; i < 5; i++) await new Promise((r) => setImmediate(r));
    expect(errors.mock.calls.map((c) => String(c[0]))).toEqual([]);
  });
});

describe('the selection marks on a trajectory frame', () => {
  async function ensemble() {
    const env = await createPlaybackPlugin(DefaultPluginSpec());
    const { plugin } = env;
    const data = await plugin.builders.data.rawData({ data: pdbEnsemble(3, 6, true), label: 'ensemble' });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
    const hierarchy = () => plugin.managers.structure.hierarchy.current;
    const frame = (i: number) => setFrame(plugin, hierarchy().models[0].cell.transform.ref, i);
    return { ...env, hierarchy, frame };
  }

  it('are not cleared and set again without a selection', async () => {
    const { canvas3d, frame } = await ensemble();
    canvas3d.mark.mockClear();
    await frame(1);
    expect(canvas3d.mark).not.toHaveBeenCalled();
  });

  it('are set again with a selection, on the new frame', async () => {
    const { plugin, canvas3d, hierarchy, frame } = await ensemble();
    const s = hierarchy().structures[0].cell.obj!.data;
    plugin.managers.interactivity.lociSelects.select({
      loci: StructureElement.Loci(s, [{ unit: s.units[0], indices: OrderedSet.ofBounds(0, 5) as any }]),
    });
    canvas3d.mark.mockClear();
    await frame(1);
    // a visual may have reset its markers (an operator change) or moved its groups (bonds, contacts of the frame): the
    // selection is marked again, on the new frame's structure
    const root = hierarchy().structures[0].cell.obj!.data;
    const selects = canvas3d.mark.mock.calls.filter((c: any[]) => c[1] === MarkerAction.Select);
    expect(selects.length).toBeGreaterThan(0);
    for (const [reprLoci] of selects) expect((reprLoci as any).loci.structure.root).toBe(root);
    const loci = plugin.managers.structure.selection.getLoci(hierarchy().structures[0].cell.obj!.data);
    expect(StructureElement.Loci.size(loci as StructureElement.Loci)).toBe(5);
  });
});

describe('a playing trajectory', () => {
  it('commits an update only on the ticks whose frame changes', async () => {
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
