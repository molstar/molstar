/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginCommands } from '../../../mol-plugin/commands';
import { StateSelection } from '../../../mol-state';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';
import { StateTransforms } from '../../transforms';
import { PluginStateAnimation } from '../model';

export const AnimateStructureSpin = PluginStateAnimation.create({
    name: 'built-in.animate-structure-spin',
    display: { name: 'Spin Structure' },
    isExportable: true,
    params: () => ({
        durationInMs: PD.Numeric(3000, { min: 100, max: 10000, step: 100 })
    }),
    initialState: () => ({ t: 0 }),
    getDuration: p => ({ kind: 'fixed', durationMs: p.durationInMs }),
    async setup(_, __, plugin) {
        const state = plugin.state.data;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D));

        const update = state.build();
        let changed = false;
        for (const r of reprs) {
            const spins = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.SpinStructureRepresentation3D, r.transform.ref));
            if (spins.length > 0) continue;

            changed = true;
            update.to(r.transform.ref)
                .apply(StateTransforms.Representation.SpinStructureRepresentation3D, { t: 0 }, { tags: 'animate-structure-spin' });
        }

        if (!changed) return;

        return update.commit({ doNotUpdateCurrent: true });
    },
    teardown(_, __, plugin) {
        const state = plugin.state.data;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3DState)
            .withTag('animate-structure-spin'));
        if (reprs.length === 0) return;

        const update = state.build();
        for (const r of reprs) update.delete(r.transform.ref);
        return update.commit();
    },
    async apply(animState, t, ctx) {
        const state = ctx.plugin.state.data;
        const anims = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.SpinStructureRepresentation3D));

        if (anims.length === 0) {
            return { kind: 'finished' };
        }

        const update = state.build();

        const d = (t.current - t.lastApplied) / ctx.params.durationInMs;
        const newTime = (animState.t + d) % 1;

        for (const m of anims) {
            update.to(m).update({ ...m.params?.values, t: newTime });
        }

        await PluginCommands.State.Update(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

        return { kind: 'next', state: { t: newTime } };
    }
});