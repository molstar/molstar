/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateObject } from '../../objects.js';
import { StateTransforms } from '../../transforms.js';
import { createTrajectoryAnimation } from '../trajectory.js';

export const AnimateModelIndex = createTrajectoryAnimation({
    name: 'built-in.animate-model-index',
    display: { name: 'Animate Trajectory' },
    transformer: StateTransforms.Model.ModelFromTrajectory,
    trajectoryType: PluginStateObject.Molecule.Trajectory,
    noTrajectoryReason: 'No trajectory to animate',
    getFrameCount: data => data.frameCount,
    getFrameIndex: params => params.modelIndex,
    setFrameIndex: modelIndex => ({ modelIndex })
});