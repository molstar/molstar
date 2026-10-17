/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { mapState } from '@molstar/prism-ui/core';
import { derive } from '@molstar/prism-ui/reactive/state';
import type { TrajectoryControlsModel } from '@molstar/prism-ui/model/playback';
import type { ViewerApi } from '../api/area.js';
import { createTrajectorySource } from '../sources/trajectory.js';

// The canvas' trajectory chip's model (prism-ui's `TrajectoryControlsModel`) over the trajectory source and the
// person's API: shown only while frames play (shell.tsx; the controls are the Structure panel's trajectory row), it
// reads the frames of the structure that plays (else the focused one's, else the first with several models) and acts
// through `trajectory.frame`, `trajectory.play` and `trajectory.stop`, never the plugin.

export function createTrajectoryControls(
  plugin: PluginContext,
  api: ViewerApi,
): { model: TrajectoryControlsModel; dispose(): void } {
  const source = createTrajectorySource(plugin);
  // The one that plays, else the one in view.
  const followed = derive([source.frames, source.playingRef, source.shown], (list, playing, shown) =>
    playing ? (list.find((t) => t.ref === playing) ?? shown) : shown,
  );
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
        if (shown) await api.trajectory.frame({ structure: shown.ref, index });
      }
    } finally {
      seeking = false;
    }
  };
  const model: TrajectoryControlsModel = {
    label: mapState(followed, (s) => s?.label ?? '', Object.is),
    count: mapState(followed, (s) => s?.count ?? 0, Object.is),
    index: mapState(followed, (s) => s?.index ?? 0, Object.is),
    isPlaying: source.playing,
    seek(index) {
      wanted = index;
      void seek();
    },
    togglePlay() {
      const shown = followed.value;
      void (source.playing.value ? api.trajectory.stop() : api.trajectory.play({ structure: shown?.ref ?? '' }));
    },
  };
  return { model, dispose: () => {} };
}
