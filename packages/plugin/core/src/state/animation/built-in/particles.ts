/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../objects.js';
import { StateTransforms } from '../../transforms.js';
import { createTrajectoryAnimation } from '../trajectory.js';

export const AnimateParticleTrajectory = createTrajectoryAnimation({
    name: 'built-in.animate-particle-trajectory',
    display: { name: 'Animate Particle Trajectory' },
    transformer: StateTransforms.Particles.ParticleListFromTrajectory,
    trajectoryType: PluginStateObject.Particle.Trajectory,
    noTrajectoryReason: 'No particle trajectory to animate',
    getFrameCount: data => data.frameCount,
    getFrameIndex: params => params.frameIndex,
    setFrameIndex: frameIndex => ({ frameIndex })
});
