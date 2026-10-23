/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { mapState, type ReadonlyState } from '@molstar/prism-ui/core';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { PeriodicBoxMode, TrajectoryControlsModel } from '@molstar/prism-ui/model/playback';
import type { ViewerApi } from '../api/area.js';
import { trajectorySourceOf } from '../sources/trajectory.js';

// The canvas' trajectory bar's model (prism-ui's `TrajectoryControlsModel`) over the shared trajectory source and the
// person's API. It follows what plays (a structure's frames, or a particle trajectory's: `particles.*`) (paused, the one that played last; else the focused one's, else the
// first with several models) and acts through `trajectory.frame`, `trajectory.play` and `trajectory.stop` (and the box
// menu through `trajectory.periodic`), never the plugin. Where it shows is the shell's (shell.tsx): a compact chip while frames play and that structure's row is on
// screen in the Structure panel, else the full bar, playing or not.

export function createTrajectoryControls(
  plugin: PluginContext,
  api: ViewerApi,
): {
  model: TrajectoryControlsModel;
  /** The structure the bar follows (its ref), and whether its Structure-panel row is on screen. */
  rowShown: ReadonlyState<boolean>;
  dispose(): void;
} {
  const source = trajectorySourceOf(plugin);
  // The one that plays; paused, the one that played last (while it is there: Pause then Play is the same one); else the
  // one in view.
  let last: string | undefined;
  const followed = derive([source.byRef, source.playingRefs, source.shown], (byRef, playing, shown) => {
    const now = playing.values().next().value;
    if (now !== undefined) last = now;
    return (last !== undefined ? byRef.get(last) : undefined) ?? shown;
  });
  // The scrubber reports every step of a drag: one frame change at a time, the last position wins.
  let seeking = false;
  let wanted: number | undefined;
  const seek = async () => {
    if (seeking) return;
    seeking = true;
    try {
      while (wanted !== undefined) {
        const index = wanted;
        wanted = undefined;
        const shown = followed.value;
        if (shown?.kind === 'particles') await api.run('particles.frame', { particles: shown.ref, index });
        else if (shown) await api.trajectory.frame({ structure: shown.ref, index });
      }
    } finally {
      seeking = false;
    }
  };
  const model: TrajectoryControlsModel = {
    label: mapState(followed, (s) => s?.label ?? '', Object.is),
    unit: mapState(followed, (s) => s?.unit ?? 'models', Object.is),
    count: mapState(followed, (s) => s?.count ?? 0, Object.is),
    index: derive([followed, source.indexByRef], (s, m) => (s ? (m.get(s.ref) ?? 0) : 0)),
    isPlaying: source.playing,
    seek(index) {
      wanted = index;
      void seek();
    },
    togglePlay() {
      const shown = followed.value;
      // A particle trajectory's frames through the particles feature (`particles.*`, its animation), else a structure's
      if (shown?.kind === 'particles') {
        void (source.playing.value ? api.run('particles.stop') : api.run('particles.play', { particles: shown.ref }));
      } else {
        void (source.playing.value ? api.trajectory.stop() : api.trajectory.play({ structure: shown?.ref ?? '' }));
      }
    },
  };
  // The box menu: the followed structure's periodic params, as one of the menu's modes
  model.periodic = {
    mode: derive([followed, source.periodicByRef], (s, byRef): PeriodicBoxMode | undefined => {
      const p = s?.kind === 'particles' ? undefined : s && byRef.get(s.ref)?.params;
      if (!p) return undefined;
      const mode = p.periodic ?? 'whole';
      return mode === 'box' && p.periodicCenter?.name === 'selection' ? 'box-centered' : mode;
    }),
    set(mode) {
      const shown = followed.value;
      if (!shown || shown.kind === 'particles') return;
      void api.trajectory.periodic(
        mode === 'box-centered'
          ? { structure: shown.ref, mode: 'box', center: 'polymers' }
          : { structure: shown.ref, mode, ...(mode === 'box' ? { center: 'box' as const } : {}) },
      );
    },
    // the neighbouring boxes through `structure.symmetry` (the cells -1 to 1 around, from the frame's box)
    images: {
      shown: derive([followed, source.periodicByRef], (s, byRef) =>
        s?.kind === 'particles' ? false : !!(s && byRef.get(s.ref)?.images),
      ),
      set(on) {
        const shown = followed.value;
        if (!shown || shown.kind === 'particles') return;
        void api.structure.symmetry(
          on ? { structure: shown.ref, kind: 'cells', cells: 1 } : { structure: shown.ref, kind: 'off' },
        );
      },
    },
  };
  const rowShown = derive([followed, source.rowsShown], (f, rows) => !!f && rows.has(f.ref));
  return { model, rowShown, dispose: () => {} };
}
