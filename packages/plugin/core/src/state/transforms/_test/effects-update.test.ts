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
import { Transparency } from '@molstar/graphics/theme/transparency';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { OverpaintStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/overpaint';
import { TransparencyStructureRepresentation3DFromScript } from '@molstar/plugin/state/transforms/structure/effects/transparency';
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
  const script = Script('(sel.atom.res (= atom.resno 2))', 'mol-script');
  const b = plugin.build();
  const overpaint = b.to(repr).apply(OverpaintStructureRepresentation3DFromScript, {
    layers: [{ script, color: 0xff0000 as any, clear: false }],
  }).ref;
  b.to(repr).apply(TransparencyStructureRepresentation3DFromScript, { layers: [{ script, value: 0.5 }] });
  await b.commit();
  return { plugin, repr, overpaint };
}

// the spies are on shared module objects
afterEach(() => vi.restoreAllMocks());

describe('an effect below an updated representation', () => {
  it('does not run its scripts again when its layers and structure are those it last ran on (S7)', async () => {
    const { plugin, repr, overpaint } = await representationWithEffects();
    const overpaintScripts = vi.spyOn(Overpaint, 'ofScript');
    const transparencyScripts = vi.spyOn(Transparency, 'ofScript');
    for (const alpha of [0.5, 0.6]) {
      await plugin
        .build()
        .to(repr)
        .update(StructureRepresentation3D, (old) => ({
          ...old,
          type: { ...old.type, params: { ...old.type.params, alpha } },
        }))
        .commit();
    }
    expect(overpaintScripts).not.toHaveBeenCalled();
    expect(transparencyScripts).not.toHaveBeenCalled();
    // the layers are still applied
    const state = plugin.state.data.cells.get(overpaint)!.obj!.data.state as { overpaint: Overpaint };
    expect(state.overpaint.layers).toHaveLength(1);
  });

  it('runs them when its layers change', async () => {
    const { plugin, overpaint } = await representationWithEffects();
    const overpaintScripts = vi.spyOn(Overpaint, 'ofScript');
    await plugin
      .build()
      .to(overpaint)
      .update(OverpaintStructureRepresentation3DFromScript, (old) => ({
        layers: [{ ...old.layers[0], color: 0x00ff00 as any }],
      }))
      .commit();
    expect(overpaintScripts).toHaveBeenCalledTimes(1);
  });
});
