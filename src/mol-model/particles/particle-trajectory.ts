/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParticleList } from './particle-list';

/**
 * A generic interface for representing (partial) particle trajectories.
 * Mirrors `Trajectory` from `mol-model/structure/trajectory.ts` but for `ParticleList` frames.
 */
export interface ParticleTrajectory {
    readonly frameCount: number

    /** Statically available representative frame. Required for example by certain UI actions. */
    readonly representative: ParticleList

    /** Return the particle list for the given frame index (0-based). */
    getFrameAtIndex(i: number): ParticleList
}

export class ArrayParticleTrajectory implements ParticleTrajectory {
    readonly frameCount: number;
    readonly representative: ParticleList;

    getFrameAtIndex(i: number): ParticleList {
        return this.frames[i];
    }

    constructor(private frames: ParticleList[]) {
        this.frameCount = frames.length;
        this.representative = frames[0];
    }
}
