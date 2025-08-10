/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginState } from '../../../mol-plugin/state';
import { PluginStateSnapshotManager } from '../../manager/snapshots';
import { PluginStateAnimation } from '../model';

async function setPartialSnapshot(plugin: PluginContext, entry: Partial<PluginStateSnapshotManager.Entry['snapshot']>, first = false) {
    if (entry.data) {
        await plugin.runTask(plugin.state.data.setSnapshot(entry.data));
        // update the canvas3d trackball with the snapshot
        plugin.canvas3d?.setProps({
            trackball: entry.canvas3d?.props?.trackball
        });

    }

    if (entry.camera?.current) {
        plugin.canvas3d?.requestCameraReset({
            snapshot: entry.camera.current,
            durationMs: first || entry.camera.transitionStyle === 'instant'
                ? 0 : entry.camera.transitionDurationInMs,
        });
    } else if (entry.camera?.focus) {
        plugin.managers.camera.focusObject({
            ...entry.camera.focus,
            durationMs: first || entry.camera.transitionStyle === 'instant'
                ? 0 : entry.camera.transitionDurationInMs,
        });
    }
}

type State = {
    totalDuration: number,
    snapshots: PluginStateSnapshotManager.Entry[],
    currentIndex: number,
    currentTransitionFrame?: number,
    isInitial?: boolean,
};

export const AnimateStateSnapshots = PluginStateAnimation.create({
    name: 'built-in.animate-state-snapshots',
    display: { name: 'State Snapshots' },
    isExportable: true,
    params: () => ({}),
    canApply(plugin) {
        const entries = plugin.managers.snapshot.state.entries;
        if (entries.size < 1) {
            return { canApply: false, reason: 'At least 1 state required.' };
        }
        if (entries.some(e => !!e?.snapshot.startAnimation)) {
            return { canApply: false, reason: 'Nested animations not supported.' };
        }
        return { canApply: plugin.managers.snapshot.state.entries.size > 0 };
    },
    setup(_, __, plugin) {
        const pivot = plugin.managers.snapshot.state.entries.get(0)!;
        setPartialSnapshot(plugin, pivot.snapshot, true);
    },
    getDuration: (_, plugin) => {
        return {
            kind: 'fixed',
            durationMs: plugin.managers.snapshot.state.entries.toArray().reduce((a, b) => a + (b.snapshot.durationInMs ?? 0), 0)
        };
    },
    initialState: (_, plugin) => {
        const snapshots = plugin.managers.snapshot.state.entries.toArray();

        return {
            totalDuration: snapshots.reduce((a, b) => a + (b.snapshot.durationInMs ?? 0), 0),
            snapshots,
            currentIndex: 0,
            currentTransitionFrame: 0,
        } as State;
    },
    async apply(animState: State, t, ctx) {
        if (t.current >= animState.totalDuration) {
            return { kind: 'finished' };
        }

        let ctime = 0, i = 0;
        let ftime = 0;
        for (const s of animState.snapshots) {
            ftime = t.current - ctime;
            ctime += s.snapshot.durationInMs ?? 0;
            if (t.current < ctime) {
                break;
            }
            i++;
        }

        if (i >= animState.snapshots.length) return { kind: 'finished' };

        const { transition, camera, canvas3d } = animState.snapshots[i].snapshot;
        const frameIndex = PluginState.getStateTransitionFrameIndex(animState.snapshots[i].snapshot, ftime);
        if (transition && frameIndex !== undefined) {
            if (i === animState.currentIndex && frameIndex === animState.currentTransitionFrame) {
                return { kind: 'skip' };
            }

            if (frameIndex === 0 || i !== animState.currentIndex) {
                await setPartialSnapshot(ctx.plugin, {
                    ...transition.frames[frameIndex],
                    camera,
                    canvas3d,
                });
            } else {
                await setPartialSnapshot(ctx.plugin, transition.frames[frameIndex]);
            }

            return { kind: 'next', state: { ...animState, currentIndex: i, currentAnimationFrame: frameIndex } };
        }

        if (i === animState.currentIndex) {
            return { kind: 'skip' };
        }

        await setPartialSnapshot(ctx.plugin, animState.snapshots[i].snapshot);
        return { kind: 'next', state: { ...animState, currentIndex: i, currentAnimationFrame: undefined } };
    }
});

export const AnimateStateSnapshotTransition = PluginStateAnimation.create({
    name: 'built-in.animate-state-snapshot-transition',
    display: { name: 'State Snapshot Transition' },
    isExportable: true,
    params: () => ({}),
    canApply(plugin) {
        const { snapshot } = plugin.managers;
        const { current } = snapshot;

        if (!current?.snapshot.transition) {
            return { canApply: false, reason: 'No transition found' };
        }

        return { canApply: true };
    },
    setup(_, __, plugin) {
        const { current } = plugin.managers.snapshot;
        if (!current) return;
        setPartialSnapshot(plugin, current.snapshot.transition?.frames[0] ?? current.snapshot, true);
    },
    getDuration: (_, plugin) => {
        const { current } = plugin.managers.snapshot;
        if (!current?.snapshot.transition) return { kind: 'fixed', durationMs: 0 };

        if (current.snapshot.transition?.loop) {
            return { kind: 'infinite' };
        }

        return {
            kind: 'fixed',
            durationMs: PluginState.getStateTransitionDuration(current.snapshot) ?? 0
        };
    },
    initialState: (_, plugin) => {
        const { current } = plugin.managers.snapshot;
        if (!current) return;

        return {
            totalDuration: current.snapshot.transition?.loop ? Number.MAX_VALUE : (PluginState.getStateTransitionDuration(current.snapshot) ?? 0),
            snapshots: [current],
            currentIndex: 0,
            currentTransitionFrame: undefined,
            isInitial: true,
        } as State;
    },
    async apply(animState: State, t, ctx) {
        const snapshot = animState.snapshots[0]?.snapshot;

        if (t.current >= animState.totalDuration) {
            if (snapshot?.transition && animState.isInitial) {
                const frameIndex = snapshot.transition.frames.length - 1;
                ctx.plugin.managers.snapshot.setSnapshotAnimationFrame(frameIndex, false);
                await setPartialSnapshot(ctx.plugin, snapshot.transition.frames[frameIndex]);
            }
            return { kind: 'finished' };
        }

        if (!snapshot) return { kind: 'finished' };

        const { transition, camera, canvas3d } = snapshot;
        const frameIndex = PluginState.getStateTransitionFrameIndex(snapshot, t.current);
        if (!transition || frameIndex === undefined) {
            return { kind: 'finished' };
        }

        if (frameIndex === animState.currentTransitionFrame) {
            return { kind: 'skip' };
        }

        ctx.plugin.managers.snapshot.setSnapshotAnimationFrame(frameIndex, false);
        if (frameIndex === 0) {
            await setPartialSnapshot(ctx.plugin, {
                ...transition.frames[frameIndex],
                camera,
                canvas3d,
            });
        } else {
            await setPartialSnapshot(ctx.plugin, transition.frames[frameIndex]);
        }
        return { kind: 'next', state: { ...animState, currentAnimationFrame: frameIndex, isInitial: false } };
    }
});