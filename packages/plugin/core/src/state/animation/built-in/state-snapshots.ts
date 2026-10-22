/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { Canvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginState } from '@molstar/plugin/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginStateSnapshotManager } from '../../manager/snapshots.js';
import { PluginStateAnimation } from '../model.js';

/** Plugins whose canvas props a playing animation set silently: its `teardown` emits `settingsUpdated` once. */
const settingsChanged = new WeakSet<PluginContext>();

/** Emits `events.canvas3d.settingsUpdated` once if frames since the last call set canvas props silently. */
function flushSettingsUpdated(plugin: PluginContext) {
  if (!settingsChanged.delete(plugin)) return;
  plugin.events.canvas3d.settingsUpdated.next(void 0);
}

/**
 * `frame`: a transition frame of a playing animation, whose canvas props are set without `settingsUpdated` (it would
 * wake every settings reader at the frame rate); the next emitting call or the animation's `teardown` emits it once.
 */
async function setPartialSnapshot(
  plugin: PluginContext,
  entry: Partial<PluginStateSnapshotManager.Entry['snapshot']>,
  first = false,
  frame = false,
) {
  if (entry.data) {
    // frames are applied at the frame rate: their timing would fill the log
    await plugin.runTask(plugin.state.data.setSnapshot(entry.data, { doNotLogTiming: true }));
  }

  if (entry.canvas3d?.props) {
    const settings = PD.normalizeParams(Canvas3DParams, entry.canvas3d.props, 'children');
    if (entry.camera?.current || entry.camera?.focus) {
      // Avoid multiple camera transitions (creates ugly cases when camera in old and new snapshot is the same)
      settings.camera = undefined;
      settings.cameraClipping = undefined;
      settings.cameraFog = undefined;
    }
    if (!frame) {
      settingsChanged.delete(plugin);
      plugin.setCanvas3DProps(settings);
    } else if (plugin.canvas3d) {
      plugin.canvas3d.setProps(settings);
      settingsChanged.add(plugin);
    }
  }
  const instant = first || entry.camera?.transitionStyle === 'instant';
  if (entry.camera?.current) {
    plugin.canvas3d?.requestCameraReset({
      snapshot: entry.camera.current,
      durationMs: instant ? 0 : entry.camera.transitionDurationInMs,
      easing: instant ? undefined : entry.camera.transitionEasing,
      trajectory: instant ? undefined : entry.camera.transitionTrajectory,
    });
  } else if (entry.camera?.focus) {
    plugin.managers.camera.focusObject({
      ...entry.camera.focus,
      durationMs: instant ? 0 : entry.camera.transitionDurationInMs,
      easing: instant ? undefined : entry.camera.transitionEasing,
      trajectory: instant ? undefined : entry.camera.transitionTrajectory,
    });
  }
}

type State = {
  totalDuration: number;
  snapshots: PluginStateSnapshotManager.Entry[];
  currentIndex: number;
  currentTransitionFrame?: number;
  isInitial?: boolean;
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
    if (entries.some((e) => !!e?.snapshot.startAnimation)) {
      return { canApply: false, reason: 'Nested animations not supported.' };
    }
    return { canApply: plugin.managers.snapshot.state.entries.size > 0 };
  },
  setup(_, __, plugin) {
    const pivot = plugin.managers.snapshot.state.entries.get(0)!;
    setPartialSnapshot(plugin, pivot.snapshot, true);
  },
  teardown(_, __, plugin) {
    flushSettingsUpdated(plugin);
  },
  getDuration: (_, plugin) => {
    return {
      kind: 'fixed',
      durationMs: plugin.managers.snapshot.state.entries
        .toArray()
        .reduce((a, b) => a + (b.snapshot.durationInMs ?? 0), 0),
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

    let ctime = 0,
      i = 0;
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
        await setPartialSnapshot(ctx.plugin, transition.frames[frameIndex], false, true);
      }

      return { kind: 'next', state: { ...animState, currentIndex: i, currentTransitionFrame: frameIndex } };
    }

    if (i === animState.currentIndex) {
      return { kind: 'skip' };
    }

    await setPartialSnapshot(ctx.plugin, animState.snapshots[i].snapshot);
    return { kind: 'next', state: { ...animState, currentIndex: i, currentTransitionFrame: undefined } };
  },
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
  teardown(_, __, plugin) {
    flushSettingsUpdated(plugin);
  },
  getDuration: (_, plugin) => {
    const { current } = plugin.managers.snapshot;
    if (!current?.snapshot.transition) return { kind: 'fixed', durationMs: 0 };

    if (current.snapshot.transition?.loop) {
      return { kind: 'infinite' };
    }

    return {
      kind: 'fixed',
      durationMs: PluginState.getStateTransitionDuration(current.snapshot) ?? 0,
    };
  },
  initialState: (_, plugin) => {
    const { current } = plugin.managers.snapshot;
    if (!current) return;

    return {
      totalDuration: current.snapshot.transition?.loop
        ? Number.MAX_VALUE
        : (PluginState.getStateTransitionDuration(current.snapshot) ?? 0),
      snapshots: [current],
      currentIndex: 0,
      currentTransitionFrame: undefined,
      isInitial: true,
    } as State;
  },
  async apply(animState: State, t, ctx) {
    const snapshot = animState.snapshots[0]?.snapshot;

    if (t.current >= animState.totalDuration) {
      // the last frame, unless it is shown already: its slot may be shorter than a tick
      const frameIndex = (snapshot?.transition?.frames.length ?? 0) - 1;
      if (snapshot?.transition && (animState.isInitial || animState.currentTransitionFrame !== frameIndex)) {
        ctx.plugin.managers.snapshot.setSnapshotAnimationFrame(animState.totalDuration, false);
        await setPartialSnapshot(ctx.plugin, snapshot.transition.frames[frameIndex], false, true);
      }
      return { kind: 'finished' };
    }

    if (!snapshot) return { kind: 'finished' };

    const { transition, camera, canvas3d } = snapshot;
    const frameIndex = PluginState.getStateTransitionFrameIndex(snapshot, t.current);
    if (!transition || frameIndex === undefined) {
      return { kind: 'finished' };
    }

    // the time moves on every tick, also while a frame holds
    ctx.plugin.managers.snapshot.setSnapshotAnimationFrame(t.current, false);
    if (frameIndex === animState.currentTransitionFrame) {
      return { kind: 'skip' };
    }

    if (frameIndex === 0) {
      await setPartialSnapshot(
        ctx.plugin,
        {
          ...transition.frames[frameIndex],
          camera,
          canvas3d,
        },
        false,
        true,
      );
    } else {
      await setPartialSnapshot(ctx.plugin, transition.frames[frameIndex], false, true);
    }
    return { kind: 'next', state: { ...animState, currentTransitionFrame: frameIndex, isInitial: false } };
  },
});
