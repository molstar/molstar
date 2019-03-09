/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateAnimation } from './model';
import { PluginStateObject } from '../objects';
import { StateTransforms } from '../transforms';
import { StateSelection } from 'mol-state';
import { PluginCommands } from 'mol-plugin/command';
import { ParamDefinition as PD } from 'mol-util/param-definition';

export const AnimateModelIndex = PluginStateAnimation.create({
    name: 'built-in.animate-model-index',
    display: { name: 'Animate Trajectory' },
    params: () => ({
        mode: PD.MappedStatic('palindrome', {
            palindrome: PD.Group({ }),
            loop: PD.Group({ }),
            once: PD.Group({ direction: PD.Select('forward', [['forward', 'Forward'], ['backward', 'Backward']]) }, { isFlat: true })
        }, { options: [['palindrome', 'Palindrome'], ['loop', 'Loop'], ['once', 'Once']] }),
        maxFPS: PD.Numeric(15, { min: 1, max: 30, step: 1 })
    }),
    initialState: () => ({} as { palindromeDirections?: { [id: string]: -1 | 1 | undefined } }),
    async apply(animState, t, ctx) {
        // limit fps
        if (t.current > 0 && t.current - t.lastApplied < 1000 / ctx.params.maxFPS) {
            return { kind: 'skip' };
        }

        const state = ctx.plugin.state.dataState;
        const models = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Model)
            .filter(c => c.transform.transformer === StateTransforms.Model.ModelFromTrajectory));

        if (models.length === 0) {
            // nothing more to do here
            return { kind: 'finished' };
        }

        const update = state.build();

        const params = ctx.params;
        const palindromeDirections = animState.palindromeDirections || { };
        let isEnd = false, allSingles = true;

        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, [PluginStateObject.Molecule.Trajectory]);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj as PluginStateObject.Molecule.Trajectory;
            update.to(m.transform.ref).update(StateTransforms.Model.ModelFromTrajectory,
                old => {
                    const len = traj.data.length;
                    if (len !== 1) {
                        allSingles = false;
                    } else {
                        return old;
                    }
                    let dir: -1 | 1 = 1;
                    if (params.mode.name === 'once') {
                        dir = params.mode.params.direction === 'backward' ? -1 : 1;
                        // if we are at start or end already, do nothing.
                        if ((dir === -1 && old.modelIndex === 0) || (dir === 1 && old.modelIndex === len - 1)) {
                            isEnd = true;
                            return old;
                        }
                    } else if (params.mode.name === 'palindrome') {
                        if (old.modelIndex === 0) dir = 1;
                        else if (old.modelIndex === len - 1) dir = -1;
                        else dir = palindromeDirections[m.transform.ref] || 1;
                    }
                    palindromeDirections[m.transform.ref] = dir;

                    let modelIndex = (old.modelIndex + dir) % len;
                    if (modelIndex < 0) modelIndex += len;

                    isEnd = isEnd || (dir === -1 && modelIndex === 0) || (dir === 1 && modelIndex === len - 1);

                    return { modelIndex };
                });
        }

        await PluginCommands.State.Update.dispatch(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

        if (allSingles || (params.mode.name === 'once' && isEnd)) return { kind: 'finished' };
        if (params.mode.name === 'palindrome') return { kind: 'next', state: { palindromeDirections } };
        return { kind: 'next', state: {} };
    }
})

export const AnimateAssemblyUnwind = PluginStateAnimation.create({
    name: 'built-in.animate-assembly-unwind',
    display: { name: 'Unwind Assembly' },
    params: () => ({
        durationInMs: PD.Numeric(3000, { min: 100, max: 10000, step: 100})
    }),
    initialState: () => ({ t: 0 }),
    async apply(animState, t, ctx) {
        const state = ctx.plugin.state.dataState;
        const anims = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure.Representation3DState)
            .filter(c => c.transform.transformer === StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D));

        if (anims.length === 0) {
            // nothing more to do here
            return { kind: 'finished' };
        }

        const update = state.build();

        const d = (t.current - t.lastApplied) / ctx.params.durationInMs;
        const newTime = (animState.t + d) % 1;

        for (const m of anims) {
            update.to(m.transform.ref).update(StateTransforms.Representation.UnwindStructureAssemblyRepresentation3D, _ => ({ t: newTime }));
        }

        await PluginCommands.State.Update.dispatch(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

        return { kind: 'next', state: { t: newTime } };
    }
})

export const AnimateUnitsExplode = PluginStateAnimation.create({
    name: 'built-in.animate-units-explode',
    display: { name: 'Explode Units' },
    params: () => ({
        durationInMs: PD.Numeric(3000, { min: 100, max: 10000, step: 100})
    }),
    initialState: () => ({ t: 0 }),
    async apply(animState, t, ctx) {
        const state = ctx.plugin.state.dataState;
        const anims = state.selectQ(q => q.rootsOfType(PluginStateObject.Molecule.Structure.Representation3DState)
            .filter(c => c.transform.transformer === StateTransforms.Representation.ExplodeStructureRepresentation3D));

        if (anims.length === 0) {
            // nothing more to do here
            return { kind: 'finished' };
        }

        const update = state.build();

        const d = (t.current - t.lastApplied) / ctx.params.durationInMs;
        const newTime = (animState.t + d) % 1;

        for (const m of anims) {
            update.to(m.transform.ref).update(StateTransforms.Representation.ExplodeStructureRepresentation3D, _ => ({ t: newTime }));
        }

        await PluginCommands.State.Update.dispatch(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });

        return { kind: 'next', state: { t: newTime } };
    }
})