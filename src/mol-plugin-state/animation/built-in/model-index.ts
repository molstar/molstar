/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginCommands } from '../../../mol-plugin/commands';
import { StateSelection } from '../../../mol-state';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';
import { StateTransforms } from '../../transforms';
import { PluginStateAnimation } from '../model';

export const AnimateModelIndex = PluginStateAnimation.create({
    name: 'built-in.animate-model-index',
    display: { name: 'Animate Trajectory' },
    params: () => ({
        mode: PD.MappedStatic('palindrome', {
            palindrome: PD.Group({ }),
            loop: PD.Group({ }),
            once: PD.Group({ direction: PD.Select('forward', [['forward', 'Forward'], ['backward', 'Backward']]) }, { isFlat: true })
        }, { options: [['palindrome', 'Palindrome'], ['loop', 'Loop'], ['once', 'Once']] }),
        maxFPS: PD.Numeric(15, { min: 1, max: 60, step: 1 })
    }),
    canApply(ctx) {
        const state = ctx.state.data;
        const models = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Model.ModelFromTrajectory));
        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, PluginStateObject.Molecule.Trajectory);
            if (parent && parent.obj && parent.obj.data.frameCount > 1) return { canApply: true };
        }
        return { canApply: false, reason: 'No trajectory to animate' };
    },
    initialState: () => ({} as { palindromeDirections?: { [id: string]: -1 | 1 | undefined } }),
    async apply(animState, t, ctx) {
        // limit fps
        if (t.current > 0 && t.current - t.lastApplied < 1000 / ctx.params.maxFPS) {
            return { kind: 'skip' };
        }

        const state = ctx.plugin.state.data;
        const models = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Model.ModelFromTrajectory));

        if (models.length === 0) {
            // nothing more to do here
            return { kind: 'finished' };
        }

        const update = state.build();

        const params = ctx.params;
        const palindromeDirections = animState.palindromeDirections || { };
        let isEnd = false, allSingles = true;

        for (const m of models) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, PluginStateObject.Molecule.Trajectory);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj;
            if (traj.data.frameCount <= 1) continue;

            update.to(m).update(old => {
                const len = traj.data.frameCount;
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

        if (!allSingles) {
            await PluginCommands.State.Update(ctx.plugin, { state, tree: update, options: { doNotLogTiming: true } });
        }

        if (allSingles || (params.mode.name === 'once' && isEnd)) return { kind: 'finished' };
        if (params.mode.name === 'palindrome') return { kind: 'next', state: { palindromeDirections } };
        return { kind: 'next', state: {} };
    }
});