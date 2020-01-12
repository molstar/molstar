/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UUID } from '../../../mol-util';
import { Cell } from '../../../mol-math/geometry/spacegroup/cell';
import { Model } from '../model';
import { AtomicConformation } from '../model/properties/atomic';
import { CustomProperties } from '../../structure';
import { Mutable } from '../../../mol-util/type-helpers';
import { Column } from '../../../mol-data/db';

export interface Frame {
    readonly elementCount: number
    readonly cell: Cell
    readonly time: Time

    // positions
    readonly x: ArrayLike<number>
    readonly y: ArrayLike<number>
    readonly z: ArrayLike<number>

    // optional velocities
    readonly velocities?: {
        readonly vx: ArrayLike<number>
        readonly vy: ArrayLike<number>
        readonly vz: ArrayLike<number>
    }

    // optional forces
    readonly forces?: {
        readonly fx: ArrayLike<number>
        readonly fy: ArrayLike<number>
        readonly fz: ArrayLike<number>
    }
}

//

export { Time }

interface Time {
    value: number
    unit: Time.Unit
}

function Time(value: number, unit: Time.Unit) {
    return { value, unit }
}

namespace Time {
    export type Unit = 'ps' | 'step'

    // TODO: conversion utilities
}

//

export { Coordinates }

interface Coordinates {
    readonly id: UUID

    readonly frames: Frame[]

    /** Number of elements (e.g. atoms) in frames */
    readonly elementCount: number

    readonly hasVelocities: boolean
    readonly hasForces: boolean

    readonly deltaTime: Time
    readonly timeOffset: Time
}

namespace Coordinates {
    export function create(frames: Frame[], deltaTime: Time, timeOffset: Time): Coordinates {
        const elementCount = frames[0].elementCount
        const hasVelocities = !!frames[0].velocities
        const hasForces = !!frames[0].forces

        return {
            id: UUID.create22(),
            frames,
            elementCount,
            hasVelocities,
            hasForces,
            deltaTime,
            timeOffset
        }
    }
}

function getAtomicConformation(frame: Frame, atomId: Column<number>): AtomicConformation {
    return {
        id: UUID.create22(),
        atomId,
        occupancy: Column.ofConst(1, frame.elementCount, Column.Schema.int),
        B_iso_or_equiv: Column.ofConst(0, frame.elementCount, Column.Schema.float),
        x: frame.x,
        y: frame.y,
        z: frame.z,
    }
}

export function trajectoryFromModelAndCoordinates(model: Model, coordinates: Coordinates): Model.Trajectory {
    const trajectory: Mutable<Model.Trajectory> = []
    const { frames } = coordinates
    for (let i = 0, il = frames.length; i < il; ++i) {
        const f = frames[i]
        const m = {
            ...model,
            id: UUID.create22(),
            modelNum: i,
            atomicConformation: getAtomicConformation(f, model.atomicConformation.atomId),
            // TODO: add support for supplying sphere and gaussian coordinates in addition to atomic coordinates
            // coarseConformation: coarse.conformation,
            customProperties: new CustomProperties(),
            _staticPropertyData: Object.create(null),
            _dynamicPropertyData: Object.create(null)
        }
        trajectory.push(m)
    }
    return trajectory
}