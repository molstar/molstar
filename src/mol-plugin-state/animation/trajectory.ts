/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginCommands } from '../../mol-plugin/commands';
import { StateObject, StateSelection, StateTransformer } from '../../mol-state';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { PluginStateAnimation } from './model';

export const TrajectoryAnimationParams = {
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
};
export type TrajectoryAnimationParams = PD.ValuesFor<typeof TrajectoryAnimationParams>

export interface TrajectoryAnimationState {
    palindromeDirections?: { [id: string]: -1 | 1 | undefined }
}

export interface TrajectoryAnimationOptions<A extends StateTransformer, TTraj extends StateObject.Ctor> {
    name: string
    display: { name: string }
    /** Transformer of the state cells that hold the current per-frame index (e.g. `ModelFromTrajectory`) */
    transformer: A
    /** Ancestor state object type that holds the trajectory data (e.g. `PluginStateObject.Molecule.Trajectory`) */
    trajectoryType: TTraj
    /** Reason shown when there is no trajectory with more than one frame to animate */
    noTrajectoryReason: string
    /** Number of frames available in the trajectory data */
    getFrameCount: (data: StateObject.From<TTraj>['data']) => number
    /** Reads the current frame index from the params of a `transformer` cell */
    getFrameIndex: (params: StateTransformer.Params<A>) => number
    /** Builds the params patch that sets a new frame index on a `transformer` cell */
    setFrameIndex: (frameIndex: number) => Partial<StateTransformer.Params<A>>
}

/**
 * Creates an animation that steps through the frames of one or more trajectories
 * (e.g. structure models or particle lists) according to the shared
 * `TrajectoryAnimationParams` (mode + duration).
 */
export function createTrajectoryAnimation<A extends StateTransformer, TTraj extends StateObject.Ctor>(options: TrajectoryAnimationOptions<A, TTraj>) {
    const { transformer, trajectoryType, noTrajectoryReason, getFrameCount, getFrameIndex, setFrameIndex } = options;

    return PluginStateAnimation.create<TrajectoryAnimationParams, TrajectoryAnimationState>({
        name: options.name,
        display: options.display,
        isExportable: true,
        params: () => TrajectoryAnimationParams,
        canApply(ctx) {
            const state = ctx.state.data;
            const cells = state.select(StateSelection.Generators.ofTransformer(transformer));
            for (const c of cells) {
                const parent = StateSelection.findAncestorOfType(state.tree, state.cells, c.transform.ref, trajectoryType);
                if (parent && parent.obj && getFrameCount(parent.obj.data) > 1) return { canApply: true };
            }
            return { canApply: false, reason: noTrajectoryReason };
        },
        getDuration: (p, ctx) => {
            if (p.duration?.name === 'fixed') {
                return { kind: 'fixed', durationMs: p.duration.params.durationInS * 1000 };
            } else if (p.duration.name === 'computed') {
                const state = ctx.state.data;
                const cells = state.select(StateSelection.Generators.ofTransformer(transformer));

                let maxDuration = 0;
                for (const c of cells) {
                    const parent = StateSelection.findAncestorOfType(state.tree, state.cells, c.transform.ref, trajectoryType);
                    if (!parent || !parent.obj) continue;
                    const frameCount = getFrameCount(parent.obj.data);
                    maxDuration = Math.max(Math.ceil(1000 * frameCount / p.duration.params.targetFps), maxDuration);
                }

                return { kind: 'fixed', durationMs: maxDuration };
            }
            return { kind: 'unknown' };
        },
        initialState: () => ({}),
        async apply(animState, t, ctx) {
            // limit fps

            if (ctx.params.duration.name === 'sequential' && t.current > 0 && t.current - t.lastApplied < 1000 / ctx.params.duration.params.maxFps) {
                return { kind: 'skip' };
            }

            const state = ctx.plugin.state.data;
            const cells = state.select(StateSelection.Generators.ofTransformer(transformer));

            if (cells.length === 0) {
                // nothing more to do here
                return { kind: 'finished' };
            }

            const update = state.build();

            const params = ctx.params;
            const palindromeDirections = animState.palindromeDirections || { };
            let isEnd = false, allSingles = true;

            for (const c of cells) {
                const parent = StateSelection.findAncestorOfType(state.tree, state.cells, c.transform.ref, trajectoryType);
                if (!parent || !parent.obj) continue;
                const len = getFrameCount(parent.obj.data);
                if (len <= 1) continue;

                update.to(c).update(old => {
                    if (len === 1) return old;
                    allSingles = false;

                    const currentIndex = getFrameIndex(old);

                    if (params.duration.name === 'sequential') {
                        let dir: -1 | 1 = 1;
                        if (params.mode.name === 'once') {
                            dir = params.mode.params.direction === 'backward' ? -1 : 1;
                            // if we are at start or end already, do nothing.
                            if ((dir === -1 && currentIndex === 0) || (dir === 1 && currentIndex === len - 1)) {
                                isEnd = true;
                                return old;
                            }
                        } else if (params.mode.name === 'palindrome') {
                            if (currentIndex === 0) dir = 1;
                            else if (currentIndex === len - 1) dir = -1;
                            else dir = palindromeDirections[c.transform.ref] || 1;
                        }
                        palindromeDirections[c.transform.ref] = dir;

                        let frameIndex = (currentIndex + dir) % len;
                        if (frameIndex < 0) frameIndex += len;

                        isEnd = isEnd || (dir === -1 && frameIndex === 0) || (dir === 1 && frameIndex === len - 1);

                        return { ...old, ...setFrameIndex(frameIndex) };
                    } else {
                        const durationInMs = params.duration.name === 'fixed'
                            ? params.duration.params.durationInS * 1000
                            : Math.ceil(1000 * len / params.duration.params.targetFps);

                        if (params.mode.name === 'once' && t.current >= durationInMs) {
                            isEnd = true;
                            return { ...old, ...setFrameIndex(len - 1) };
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

                        const frameIndex = Math.min(Math.floor(len * phase), len - 1);
                        return { ...old, ...setFrameIndex(frameIndex) };
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
}
