/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '@molstar/plugin/commands';
import { StateSelection } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginStateObject } from '../../objects.js';
import { ExplodeStructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/animation';
import { PluginStateAnimation } from '../model.js';

export const AnimateUnitsExplode = PluginStateAnimation.create({
  name: 'built-in.animate-units-explode',
  display: { name: 'Explode Units' },
  params: () => ({
    durationInMs: PD.Numeric(3000, { min: 100, max: 10000, step: 100 }),
  }),
  initialState: () => ({ t: 0 }),
  async setup(_, __, plugin) {
    const state = plugin.state.data;
    const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D));

    const update = state.build();
    let changed = false;
    for (const r of reprs) {
      const explodes = state.select(
        StateSelection.Generators.ofTransformer(ExplodeStructureRepresentation3D, r.transform.ref),
      );
      if (explodes.length > 0) continue;

      changed = true;
      update.to(r.transform.ref).apply(ExplodeStructureRepresentation3D, { t: 0 }, { tags: 'animate-units-explode' });
    }

    if (!changed) return;

    return update.commit({ doNotUpdateCurrent: true });
  },
  teardown(_, __, plugin) {
    const state = plugin.state.data;
    const reprs = state.select(
      StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3DState).withTag(
        'animate-units-explode',
      ),
    );
    if (reprs.length === 0) return;

    const update = state.build();
    for (const r of reprs) update.delete(r.transform.ref);
    return update.commit();
  },
  async apply(animState, t, ctx) {
    const state = ctx.plugin.state.data;
    const anims = state.select(StateSelection.Generators.ofTransformer(ExplodeStructureRepresentation3D));

    if (anims.length === 0) {
      return { kind: 'finished' };
    }

    const update = state.build();

    const d = (t.current - t.lastApplied) / ctx.params.durationInMs;
    const newTime = (animState.t + d) % 1;

    for (const m of anims) {
      update.to(m).update({ t: newTime });
    }

    await PluginCommands.State.Update(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

    return { kind: 'next', state: { t: newTime } };
  },
});
