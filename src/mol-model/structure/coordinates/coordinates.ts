/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UUID } from '../../../mol-util';
import { Cell } from '../../../mol-math/geometry/spacegroup/cell';
import { AtomicConformation } from '../model/properties/atomic';
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

export { Time };

interface Time {
    value: number
    unit: Time.Unit
}

function Time(value: number, unit: Time.Unit) {
    return { value, unit };
}

namespace Time {
    export type Unit = 'ps' | 'step'

    // TODO: conversion utilities
}

//

export { Coordinates };

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
        const elementCount = frames[0].elementCount;
        const hasVelocities = !!frames[0].velocities;
        const hasForces = !!frames[0].forces;

        return {
            id: UUID.create22(),
            frames,
            elementCount,
            hasVelocities,
            hasForces,
            deltaTime,
            timeOffset
        };
    }

    export function getAtomicConformation(frame: Frame, atomId: Column<number>): AtomicConformation {
        return {
            id: UUID.create22(),
            atomId,
            occupancy: Column.ofConst(1, frame.elementCount, Column.Schema.int),
            B_iso_or_equiv: Column.ofConst(0, frame.elementCount, Column.Schema.float),
            xyzDefined: true,
            x: frame.x,
            y: frame.y,
            z: frame.z,
        };
    }

    function reorderCoords(xs: ArrayLike<number>, index: ArrayLike<number>) {
        const ret = new Float32Array(xs.length);
        for (let i = 0, _i = xs.length; i < _i; i++) {
            ret[i] = xs[index[i]];
        }
        return ret;
    }

    export function getAtomicConformationReordered(frame: Frame, atomId: Column<number>, srcIndex: ArrayLike<number>): AtomicConformation {
        return {
            id: UUID.create22(),
            atomId,
            occupancy: Column.ofConst(1, frame.elementCount, Column.Schema.int),
            B_iso_or_equiv: Column.ofConst(0, frame.elementCount, Column.Schema.float),
            xyzDefined: true,
            x: reorderCoords(frame.x, srcIndex),
            y: reorderCoords(frame.y, srcIndex),
            z: reorderCoords(frame.z, srcIndex)
        };
    }
}