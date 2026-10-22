/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// An effect cell (overpaint, transparency, ...) below a representation is evaluated again whenever the representation
// is (its child): a snapshot frame that changes the representation's own params, for example.

import { afterEach, describe, expect, it, vi } from 'vitest';
import { Script } from '@molstar/model/script/script';
import { Overpaint } from '@molstar/graphics/theme/overpaint';
import { OrderedSet } from '@molstar/core/data/int';
import { StructureElement } from '@molstar/model/model/structure';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import {
  OverpaintStructureRepresentation3DFromBundle,
  OverpaintStructureRepresentation3DFromScript,
} from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { pdbEnsemble } from '../../animation/_test/harness.js';

async function representationWithEffects() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  const data = await plugin.builders.data.rawData({ data: pdbEnsemble(1, 6), label: 'chain' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  const repr = await plugin.builders.structure.representation.addRepresentation(structure, { type: 'ball-and-stick' });
  const s = structure.obj!.data;
  const bundle = StructureElement.Bundle.fromLoci(
    StructureElement.Loci(s, [{ unit: s.units[0], indices: OrderedSet.ofBounds(5, 10) as any }]),
  );
  const script = Script('(sel.atom.res (= atom.resno 2))', 'mol-script');
  const b = plugin.build();
  const fromBundle = b.to(repr).apply(OverpaintStructureRepresentation3DFromBundle, {
    layers: [{ bundle, color: 0xff0000 as any, clear: false }],
  }).ref;
  const fromScript = b.to(repr).apply(OverpaintStructureRepresentation3DFromScript, {
    layers: [{ script, color: 0x00ff00 as any, clear: false }],
  }).ref;
  await b.commit();
  const updateAlpha = (alpha: number) =>
    plugin
      .build()
      .to(repr)
      .update(StructureRepresentation3D, (old) => ({
        ...old,
        type: { ...old.type, params: { ...old.type.params, alpha } },
      }))
      .commit();
  return { plugin, fromBundle, fromScript, updateAlpha };
}

// the spies are on shared module objects
afterEach(() => vi.restoreAllMocks());

describe('an effect below an updated representation', () => {
  it('does not apply its bundles again when they and the structure are those it last ran on (S7)', async () => {
    const { plugin, fromBundle, updateAlpha } = await representationWithEffects();
    const bundles = vi.spyOn(Overpaint, 'ofBundle');
    for (const alpha of [0.5, 0.6]) await updateAlpha(alpha);
    expect(bundles).not.toHaveBeenCalled();
    // the layers are still applied
    const state = plugin.state.data.cells.get(fromBundle)!.obj!.data.state as { overpaint: Overpaint };
    expect(state.overpaint.layers).toHaveLength(1);
  });

  it('runs its scripts again: a script may read a property an ancestor computed again on the same structure', async () => {
    const { updateAlpha } = await representationWithEffects();
    const scripts = vi.spyOn(Overpaint, 'ofScript');
    for (const alpha of [0.5, 0.6]) await updateAlpha(alpha);
    expect(scripts).toHaveBeenCalledTimes(2);
  });

  it('applies its layers when they change', async () => {
    const { plugin, fromBundle } = await representationWithEffects();
    const bundles = vi.spyOn(Overpaint, 'ofBundle');
    await plugin
      .build()
      .to(fromBundle)
      .update(OverpaintStructureRepresentation3DFromBundle, (old) => ({
        layers: [{ ...old.layers[0], color: 0x0000ff as any }],
      }))
      .commit();
    expect(bundles).toHaveBeenCalledTimes(1);
  });
});
