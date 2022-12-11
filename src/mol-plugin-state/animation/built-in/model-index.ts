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
    isExportable: true,
    params: () => ({
        mode: PD.MappedStatic('loop', {
            palindrome: PD.Group({ }),
            loop: PD.Group({ direction: PD.Select('forward', [['forward', 'Forward'], ['backward', 'Backward']]) }),
            once: PD.Group({ direction: PD.Select('forward', [['forward', 'Forward'], ['backward', 'Backward']]) }, { isFlat: true })
        }, { options: [['palindrome', 'Palindrome'], ['loop', 'Loop'], ['once', 'Once']] }),
        duration: PD.MappedStatic('fixed', {
            fixed: PD.Group({
                durationInS: PD.Numeric(5, { min: 1, max: 120, step: 0.1 }, { description: 'Duration in seconds' })
            }, { isFlat: true }),
            computed: PD.Group({
                targetFps: PD.Numeric(30, { min: 5, max: 250, step: 1 }, { label: 'Target FPS' })
            }, { isFlat: true }),
            sequential: PD.Group({
                maxFps: PD.Numeric(30, { min: 5, max: 60, step: 1 })
            }, { isFlat: true })
        })
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
    getDuration: (p, ctx) => {
        if (p.duration?.name === 'fixed') {
            return { kind: 'fixed', durationMs: p.duration.params.durationInS * 1000 };
        } else if (p.duration.name === 'computed') {
            const state = ctx.state.data;
            const models = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Model.ModelFromTrajectory));

            let maxDuration = 0;
            for (const m of models) {
                const parent = StateSelection.findAncestorOfType(state.tree, state.cells, m.transform.ref, PluginStateObject.Molecule.Trajectory);
                if (!parent || !parent.obj) continue;
                const traj = parent.obj;
                maxDuration = Math.max(Math.ceil(1000 * traj.data.frameCount / p.duration.params.targetFps), maxDuration);
            }

            return { kind: 'fixed', durationMs: maxDuration };
        }
        return { kind: 'unknown' };
    },
    initialState: () => ({} as { palindromeDirections?: { [id: string]: -1 | 1 | undefined } }),
    async apply(animState, t, ctx) {
        // limit fps

        if (ctx.params.duration.name === 'sequential' && t.current > 0 && t.current - t.lastApplied < 1000 / ctx.params.duration.params.maxFps) {
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

                if (params.duration.name === 'sequential') {
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
                } else {
                    const durationInMs = params.duration.name === 'fixed'
                        ? params.duration.params.durationInS * 1000
                        : Math.ceil(1000 * traj.data.frameCount / params.duration.params.targetFps);

                    if (params.mode.name === 'once' && t.current >= durationInMs) {
                        isEnd = true;
                        return { modelIndex: traj.data.frameCount - 1 };
                    }

                    let phase: number = (t.current % durationInMs) / durationInMs;
                    if (params.mode.name === 'loop') {
                        if (params.mode.params.direction === 'backward') {
                            phase = 1 - phase;
                        }
                    } if (params.mode.name === 'palindrome') {
                        phase = 2 * phase;
                        if (phase > 1) phase = 2 - phase;
                    }

                    const modelIndex = Math.min(Math.floor(traj.data.frameCount * phase), traj.data.frameCount - 1);
                    return { modelIndex };
                }
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