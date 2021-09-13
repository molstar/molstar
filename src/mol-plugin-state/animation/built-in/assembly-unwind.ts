/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateAnimation } from '../model';
import { PluginStateObject } from '../../objects';
import { StateTransforms } from '../../transforms';
import { StateSelection, StateTransform } from '../../../mol-state';
import { PluginCommands } from '../../../mol-plugin/commands';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginContext } from '../../../mol-plugin/context';

export const AnimateAssemblyUnwind = PluginStateAnimation.create({
    name: 'built-in.animate-assembly-unwind',
    display: { name: 'Unwind Assembly' },
    isExportable: true,
    params: (plugin: PluginContext) => {
        const targets: [string, string][] = [['all', 'All']];
        const structures = plugin.state.data.select(StateSelection.Generators.rootsOfType(PluginStateObject.Molecule.Structure));

        for (const s of structures) {
            targets.push([s.transform.ref, s.obj!.data.models[0].label]);
        }

        return {
            durationInMs: PD.Numeric(3000, { min: 100, max: 10000, step: 100 }),
            playOnce: PD.Boolean(false),
            target: PD.Select(targets[0][0], targets)
        };
    },
    canApply(plugin) {
        const state = plugin.state.data;
        const root = StateTransform.RootRef;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D, root));
        return { canApply: reprs.length > 0 };
    },
    getDuration: (params) => {
        return {
            kind: 'fixed',
            durationMs: params.durationInMs
        };
    },
    initialState: () => ({ t: 0 }),
    setup(params, _, plugin) {
        const state = plugin.state.data;
        const root = !params.target || params.target === 'all' ? StateTransform.RootRef : params.target;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3D, root));

        const update = state.build();
        let changed = false;
        for (const r of reprs) {
            const unwinds = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D, r.transform.ref));
            if (unwinds.length > 0) continue;

            changed = true;
            update.to(r)
                .apply(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D, { t: 0 }, { tags: 'animate-assembly-unwind' });
        }

        if (!changed) return;

        return update.commit({ doNotUpdateCurrent: true });
    },
    teardown(_, __, plugin) {
        const state = plugin.state.data;
        const reprs = state.select(StateSelection.Generators.ofType(PluginStateObject.Molecule.Structure.Representation3DState)
            .withTag('animate-assembly-unwind'));
        if (reprs.length === 0) return;

        const update = state.build();
        for (const r of reprs) update.delete(r.transform.ref);
        return update.commit();
    },
    async apply(animState, t, ctx) {
        const state = ctx.plugin.state.data;
        const root = !ctx.params.target || ctx.params.target === 'all' ? StateTransform.RootRef : ctx.params.target;
        const anims = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D, root));

        if (anims.length === 0) {
            return { kind: 'finished' };
        }

        const update = state.build();

        const d = (t.current - t.lastApplied) / ctx.params.durationInMs;
        let newTime = (animState.t + d), finished = false;
        if (ctx.params.playOnce && newTime >= 1) {
            finished = true;
            newTime = 1;
        } else {
            newTime = newTime % 1;
        }

        for (const m of anims) {
            update.to(m).update({ t: newTime });
        }

        await PluginCommands.State.Update(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

        if (finished) return { kind: 'finished' };
        return { kind: 'next', state: { t: newTime } };
    }
});