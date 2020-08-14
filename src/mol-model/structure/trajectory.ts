/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../../mol-task';
import { Model } from '../structure';

export type TrajectoryFrameType =
  | { type: 'default' }
  /** Returns the closest available frame to the requested index  */
  | { type: 'snap' }
  /** Interpolates between two available adjacent frames */
  | { type: 'interpolate', kind?: 'linear' }

/**
 * A generic interface for representing (partial) trajectories
 */
export interface Trajectory {
    readonly duration: number,
    readonly frameCount: number,

    /** Statically available representative model. Required for example by certain UI actions. */
    readonly representative: Model,

    /** Allows to asynchronously query data from a server or interpolate frames on the fly */
    getFrameAtIndex(i: number, type?: TrajectoryFrameType): Task<Model> | Model
}

export class ArrayTrajectory implements Trajectory {
    readonly duration: number;
    readonly frameCount: number;
    readonly representative: Model;

    getFrameAtIndex(i: number) {
        return this.frames[i];
    }

    constructor(private frames: Model[]) {
        this.frameCount = frames.length;
        this.representative = frames[0];
        this.duration = frames.length;
    }
}