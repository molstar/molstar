/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { TrajectoryApi } from '../../api/trajectory.js';
import { trajectoryCommands } from '../../commands/trajectory.js';
import type { ViewerExtension } from '../../spec.js';

// The `trajectory` feature (viewer-spec.md §3.1): a trajectory's frames, played and stepped (`trajectory`), and
// its commands; the shell's trajectory bar shows while it is here and a structure has frames.

export const TrajectoryFeature: ViewerExtension = {
  id: 'trajectory',
  title: 'Trajectory',
  builtIn: true,
  api: { trajectory: (c) => new TrajectoryApi(c) },
  commands: (api) => trajectoryCommands(api),
};
