/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayEqual, UUID } from '../../../mol-util';
import { Cell } from '../../../mol-math/geometry/spacegroup/cell';
import { AtomicConformation } from '../model/properties/atomic';
import { Column } from '../../../mol-data/db';

export interface Frame {
    readonly elementCount: number
    readonly time: Time

    // positions
    readonly x: ArrayLike<number>
    readonly y: ArrayLike<number>
    readonly z: ArrayLike<number>

    // optional cell
    readonly cell?: Cell

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

    readonly xyzOrdering: {
        isIdentity: boolean,
        frozen?: boolean,
        index?: ArrayLike<number>,
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

    readonly hasCell: boolean
    readonly hasVelocities: boolean
    readonly hasForces: boolean

    readonly deltaTime: Time
    readonly timeOffset: Time
}

namespace Coordinates {
    export function create(frames: Frame[], deltaTime: Time, timeOffset: Time): Coordinates {
        const hasCell = !!frames[0].cell;
        const hasVelocities = !!frames[0].velocities;
        const hasForces = !!frames[0].forces;

        return {
            id: UUID.create22(),
            frames,
            hasCell,
            hasVelocities,
            hasForces,
            deltaTime,
            timeOffset,
        };
    }

    /**
     * Only use ordering if it's not identity.
     */
    export function getAtomicConformation(frame: Frame, fields: { atomId: Column<number>, occupancy?: Column<number>, B_iso_or_equiv?: Column<number> }, ordering?: ArrayLike<number>): AtomicConformation {
        let { x, y, z } = frame;

        if (frame.xyzOrdering.frozen) {
            if (ordering) {
                if (frame.xyzOrdering.isIdentity) {
                    // simple list reordering
                    x = getOrderedCoords(x, ordering);
                    y = getOrderedCoords(y, ordering);
                    z = getOrderedCoords(z, ordering);
                } else if (!arrayEqual(frame.xyzOrdering.index! as any, ordering as any)) {
                    x = getSourceOrderedCoords(x, frame.xyzOrdering.index!, ordering);
                    y = getSourceOrderedCoords(y, frame.xyzOrdering.index!, ordering);
                    z = getSourceOrderedCoords(z, frame.xyzOrdering.index!, ordering);
                }
            } else if (!frame.xyzOrdering.isIdentity) {
                x = getInvertedCoords(x, frame.xyzOrdering.index!);
                y = getInvertedCoords(y, frame.xyzOrdering.index!);
                z = getInvertedCoords(z, frame.xyzOrdering.index!);
            }
        } else if (ordering) {
            if (frame.xyzOrdering.isIdentity) {
                frame.xyzOrdering.isIdentity = false;
                frame.xyzOrdering.index = ordering;
                reorderCoordsInPlace(x as unknown as number[], ordering);
                reorderCoordsInPlace(y as unknown as number[], ordering);
                reorderCoordsInPlace(z as unknown as number[], ordering);
            } else {
                // is current ordering is not the same as requested?
                //   => copy the conformations into a new array
                if (!arrayEqual(frame.xyzOrdering.index! as any, ordering as any)) {
                    x = getSourceOrderedCoords(x, frame.xyzOrdering.index!, ordering);
                    y = getSourceOrderedCoords(y, frame.xyzOrdering.index!, ordering);
                    z = getSourceOrderedCoords(z, frame.xyzOrdering.index!, ordering);
                }
            }
        }

        // once the conformation has been accessed at least once, freeze it.
        //   => any other request to the frame with different ordering will result in a copy.
        frame.xyzOrdering.frozen = true;

        return {
            id: UUID.create22(),
            atomId: fields.atomId,
            occupancy: fields.occupancy ?? Column.ofConst(1, frame.elementCount, Column.Schema.int),
            B_iso_or_equiv: fields.B_iso_or_equiv ?? Column.ofConst(0, frame.elementCount, Column.Schema.float),
            xyzDefined: true,
            x,
            y,
            z,
        };
    }

    const _reorderBuffer = [0.123];
    function reorderCoordsInPlace(xs: number[], index: ArrayLike<number>) {
        const buffer = _reorderBuffer;

        for (let i = 0, _i = xs.length; i < _i; i++) {
            buffer[i] = xs[index[i]];
        }
        for (let i = 0, _i = xs.length; i < _i; i++) {
            xs[i] = buffer[i];
        }
    }

    function getSourceOrderedCoords(xs: ArrayLike<number>, srcIndex: ArrayLike<number>, index: ArrayLike<number>) {
        const ret = new Float32Array(xs.length);

        for (let i = 0, _i = xs.length; i < _i; i++) {
            ret[i] = xs[srcIndex[index[i]]];
        }

        return ret;
    }

    function getOrderedCoords(xs: ArrayLike<number>, index: ArrayLike<number>) {
        const ret = new Float32Array(xs.length);

        for (let i = 0, _i = xs.length; i < _i; i++) {
            ret[i] = xs[index[i]];
        }

        return ret;
    }

    function getInvertedCoords(xs: ArrayLike<number>, index: ArrayLike<number>) {
        const ret = new Float32Array(xs.length);

        for (let i = 0, _i = xs.length; i < _i; i++) {
            ret[index[i]] = xs[i];
        }

        return ret;
    }
}