/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { PluginStateObject } from '../../objects';
import { StateTransforms } from '../../transforms';
import { createTrajectoryAnimation } from '../trajectory';

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
