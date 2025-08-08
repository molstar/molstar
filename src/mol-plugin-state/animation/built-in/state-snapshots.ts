/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '../../../mol-plugin/context';
import { PluginStateSnapshotManager } from '../../manager/snapshots';
import { PluginStateAnimation } from '../model';

async function setPartialSnapshot(plugin: PluginContext, snapshot: Partial<PluginStateSnapshotManager.Entry['snapshot']>, first = false) {
    if (snapshot.data) {
        console.time('snapshot')
        await plugin.runTask(plugin.state.data.setSnapshot(snapshot.data));
        console.timeEnd('snapshot')
        // update the canvas3d trackball with the snapshot
        plugin.canvas3d?.setProps({
            trackball: snapshot.canvas3d?.props?.trackball
        });

    }

    // if (snapshot.camera?.current) {
    //     plugin.canvas3d?.requestCameraReset({
    //         snapshot: snapshot.camera.current,
    //         durationMs: first || snapshot.camera.transitionStyle === 'instant'
    //             ? 0 : snapshot.camera.transitionDurationInMs,
    //     });
    // } else if (snapshot.camera?.focus) {
    //     plugin.managers.camera.focusObject({
    //         ...snapshot.camera.focus,
    //         durationMs: first || snapshot.camera.transitionStyle === 'instant'
    //             ? 0 : snapshot.camera.transitionDurationInMs,
    //     });
    // }
}

type State = { totalDuration: number, snapshots: PluginStateSnapshotManager.Entry[], currentIndex: number, transitionIndex: number };

export const AnimateStateSnapshots = PluginStateAnimation.create({
    name: 'built-in.animate-state-snapshots',
    display: { name: 'State Snapshots' },
    isExportable: true,
    params: () => ({}),
    canApply(plugin) {
        const entries = plugin.managers.snapshot.state.entries;
        if (entries.size < 2) {
            return { canApply: false, reason: 'At least 2 states required.' };
        }
        if (entries.some(e => !!e?.snapshot.startAnimation)) {
            return { canApply: false, reason: 'Nested animations not supported.' };
        }
        return { canApply: plugin.managers.snapshot.state.entries.size > 1 };
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
            transitionIndex: 0,
        } as State;
    },
    async apply(animState: State, t, ctx) {
        if (t.current >= animState.totalDuration) {
            return { kind: 'finished' };
        }

        let ctime = 0, i = 0, ptime = 0;
        for (const s of animState.snapshots) {
            ctime += s.snapshot.durationInMs ?? 0;
            if (t.current < ctime) {
                break;
            }
            ptime += s.snapshot.durationInMs ?? 0;
            i++;
        }

        if (i >= animState.snapshots.length) return { kind: 'finished' };

        const transition = animState.snapshots[i].snapshot.transition;
        if (transition) {
            const stime = t.current - ptime;
            const dt = transition.durationInMs / ((transition.frames.length) || 1);
            const frame = Math.min(Math.floor(stime / dt), transition.frames.length - 1);
            if (i === animState.currentIndex && frame === animState.transitionIndex) {
                return { kind: 'skip' };
            }
            // TODO: skip if current
            await setPartialSnapshot(ctx.plugin, transition.frames[frame]);
            return { kind: 'next', state: { ...animState, currentIndex: i } };
        }

        if (i === animState.currentIndex) {
            return { kind: 'skip' };
        }
        await setPartialSnapshot(ctx.plugin, animState.snapshots[i].snapshot);

        return { kind: 'next', state: { ...animState, currentIndex: i, transitionIndex: 0 } };
    }
});