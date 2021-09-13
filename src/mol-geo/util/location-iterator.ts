/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Iterator } from '../../mol-data';
import { Vec3 } from '../../mol-math/linear-algebra';
import { NullLocation, Location } from '../../mol-model/location';

export interface LocationValue {
    location: Location
    index: number
    groupIndex: number
    instanceIndex: number
    isSecondary: boolean
}

export interface LocationIterator extends Iterator<LocationValue> {
    readonly hasNext: boolean
    readonly isNextNewInstance: boolean
    readonly groupCount: number
    readonly instanceCount: number
    readonly count: number
    readonly stride: number
    readonly nonInstanceable: boolean
    move(): LocationValue
    reset(): void
    skipInstance(): void
    voidInstances(): void
}

type LocationGetter = (groupIndex: number, instanceIndex: number) => Location
type IsSecondaryGetter = (groupIndex: number, instanceIndex: number) => boolean

export function LocationIterator(groupCount: number, instanceCount: number, stride: number, getLocation: LocationGetter, nonInstanceable = false, isSecondary: IsSecondaryGetter = () => false): LocationIterator {
    if (groupCount % stride !== 0) {
        throw new Error('incompatible groupCount and stride');
    }

    const value: LocationValue = {
        location: NullLocation as Location,
        index: 0,
        groupIndex: 0,
        instanceIndex: 0,
        isSecondary: false
    };

    let hasNext = value.groupIndex < groupCount;
    let isNextNewInstance = false;
    let groupIndex = 0;
    let instanceIndex = 0;
    let voidInstances = false;

    return {
        get hasNext() { return hasNext; },
        get isNextNewInstance() { return isNextNewInstance; },
        groupCount,
        instanceCount,
        count: groupCount * instanceCount,
        stride,
        nonInstanceable,
        move() {
            if (hasNext) {
                value.groupIndex = groupIndex;
                value.instanceIndex = instanceIndex;
                value.index = instanceIndex * groupCount + groupIndex;
                value.location = getLocation(groupIndex, voidInstances ? -1 : instanceIndex);
                value.isSecondary = isSecondary(groupIndex, voidInstances ? -1 : instanceIndex);
                groupIndex += stride;
                if (groupIndex === groupCount) {
                    ++instanceIndex;
                    isNextNewInstance = true;
                    if (instanceIndex < instanceCount) groupIndex = 0;
                } else {
                    isNextNewInstance = false;
                }
                hasNext = groupIndex < groupCount;
            }
            return value;
        },
        reset() {
            value.location = NullLocation;
            value.index = 0;
            value.groupIndex = 0;
            value.instanceIndex = 0;
            value.isSecondary = false;

            hasNext = value.groupIndex < groupCount;
            isNextNewInstance = false;
            groupIndex = 0;
            instanceIndex = 0;
            voidInstances = false;
        },
        skipInstance() {
            if (hasNext && value.instanceIndex === instanceIndex) {
                ++instanceIndex;
                groupIndex = 0;
                hasNext = instanceIndex < instanceCount;
            }
        },
        voidInstances() {
            voidInstances = true;
        }
    };
}

//

/** A position Location */
export interface PositionLocation {
    readonly kind: 'position-location',
    readonly position: Vec3
}
export function PositionLocation(position?: Vec3): PositionLocation {
    return { kind: 'position-location', position: position ? Vec3.clone(position) : Vec3() };
}
export function isPositionLocation(x: any): x is PositionLocation {
    return !!x && x.kind === 'position-location';
}