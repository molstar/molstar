/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginContext } from '@molstar/plugin/context';
import { mapState } from '@molstar/prism-ui/core';
import type { TrajectoryControlsModel } from '@molstar/prism-ui/model/playback';
import type { ViewerApi } from '../api/area.js';
import { createTrajectorySource } from '../sources/trajectory.js';

// The trajectory bar's model (prism-ui's `TrajectoryControlsModel`) over the trajectory source and the person's API:
// it reads the frames of the structure in view (the focused one's, else the first with several models) and acts
// through `trajectory.frame`, `trajectory.play` and `trajectory.stop`, never the plugin.

export function createTrajectoryControls(
  plugin: PluginContext,
  api: ViewerApi,
): { model: TrajectoryControlsModel; dispose(): void } {
  const source = createTrajectorySource(plugin);
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
        const shown = source.shown.value;
        if (shown) await api.trajectory.frame({ structure: shown.ref, index });
      }
    } finally {
      seeking = false;
    }
  };
  const model: TrajectoryControlsModel = {
    label: mapState(source.shown, (s) => s?.label ?? '', Object.is),
    count: mapState(source.shown, (s) => s?.count ?? 0, Object.is),
    index: mapState(source.shown, (s) => s?.index ?? 0, Object.is),
    isPlaying: source.playing,
    seek(index) {
      wanted = index;
      void seek();
    },
    togglePlay() {
      void (source.playing.value ? api.trajectory.stop() : api.trajectory.play());
    },
  };
  return { model, dispose: () => {} };
}
