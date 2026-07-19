/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginCommands } from '../../../mol-plugin/commands';
import { StateSelection } from '../../../mol-state';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { PluginStateObject } from '../../objects';
import { StateTransforms } from '../../transforms';
import { PluginStateAnimation } from '../model';

export const AnimateParticleTrajectory = PluginStateAnimation.create({
    name: 'built-in.animate-particle-trajectory',
    display: { name: 'Animate Particle Trajectory' },
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
        const lists = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Particles.ParticleListFromTrajectory));
        for (const l of lists) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, l.transform.ref, PluginStateObject.Particle.Trajectory);
            if (parent && parent.obj && parent.obj.data.frameCount > 1) return { canApply: true };
        }
        return { canApply: false, reason: 'No particle trajectory to animate' };
    },
    getDuration: (p, ctx) => {
        if (p.duration?.name === 'fixed') {
            return { kind: 'fixed', durationMs: p.duration.params.durationInS * 1000 };
        } else if (p.duration.name === 'computed') {
            const state = ctx.state.data;
            const lists = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Particles.ParticleListFromTrajectory));

            let maxDuration = 0;
            for (const l of lists) {
                const parent = StateSelection.findAncestorOfType(state.tree, state.cells, l.transform.ref, PluginStateObject.Particle.Trajectory);
                if (!parent || !parent.obj) continue;
                maxDuration = Math.max(Math.ceil(1000 * parent.obj.data.frameCount / p.duration.params.targetFps), maxDuration);
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
        const lists = state.select(StateSelection.Generators.ofTransformer(StateTransforms.Particles.ParticleListFromTrajectory));

        if (lists.length === 0) {
            return { kind: 'finished' };
        }

        const update = state.build();

        const params = ctx.params;
        const palindromeDirections = animState.palindromeDirections || { };
        let isEnd = false, allSingles = true;

        for (const l of lists) {
            const parent = StateSelection.findAncestorOfType(state.tree, state.cells, l.transform.ref, PluginStateObject.Particle.Trajectory);
            if (!parent || !parent.obj) continue;
            const traj = parent.obj;
            if (traj.data.frameCount <= 1) continue;

            update.to(l).update(old => {
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
                        if ((dir === -1 && old.frameIndex === 0) || (dir === 1 && old.frameIndex === len - 1)) {
                            isEnd = true;
                            return old;
                        }
                    } else if (params.mode.name === 'palindrome') {
                        if (old.frameIndex === 0) dir = 1;
                        else if (old.frameIndex === len - 1) dir = -1;
                        else dir = palindromeDirections[l.transform.ref] || 1;
                    }
                    palindromeDirections[l.transform.ref] = dir;

                    let frameIndex = (old.frameIndex + dir) % len;
                    if (frameIndex < 0) frameIndex += len;

                    isEnd = isEnd || (dir === -1 && frameIndex === 0) || (dir === 1 && frameIndex === len - 1);

                    return { frameIndex };
                } else {
                    const durationInMs = params.duration.name === 'fixed'
                        ? params.duration.params.durationInS * 1000
                        : Math.ceil(1000 * traj.data.frameCount / params.duration.params.targetFps);

                    if (params.mode.name === 'once' && t.current >= durationInMs) {
                        isEnd = true;
                        return { frameIndex: traj.data.frameCount - 1 };
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

                    const frameIndex = Math.min(Math.floor(traj.data.frameCount * phase), traj.data.frameCount - 1);
                    return { frameIndex };
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
