/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Cai Huiyu <szmun.caihy@gmail.com>
 */

import { Iterator } from '../../mol-data';
import { Vec3 } from '../../mol-math/linear-algebra';
import { NullLocation, Location } from '../../mol-model/location';

export interface LocationValue {
    location: Location
    location2: Location
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
    readonly hasLocation2: boolean
    move(): LocationValue
    reset(): void
    skipInstance(): void
    voidInstances(): void
}

type LocationGetter = (groupIndex: number, instanceIndex: number) => Location
type IsSecondaryGetter = (groupIndex: number, instanceIndex: number) => boolean

export function LocationIterator(groupCount: number, instanceCount: number, stride: number, getLocation: LocationGetter, nonInstanceable = false, isSecondary: IsSecondaryGetter = () => false, getLocation2?: LocationGetter): LocationIterator {
    if (groupCount % stride !== 0) {
        throw new Error('incompatible groupCount and stride');
    }

    const value: LocationValue = {
        location: NullLocation as Location,
        location2: NullLocation as Location,
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

    const hasLocation2 = !!getLocation2;

    return {
        get hasNext() { return hasNext; },
        get isNextNewInstance() { return isNextNewInstance; },
        groupCount,
        instanceCount,
        count: groupCount * instanceCount,
        stride,
        nonInstanceable,
        hasLocation2,
        move() {
            if (hasNext) {
                value.groupIndex = groupIndex;
                value.instanceIndex = instanceIndex;
                value.index = instanceIndex * groupCount + groupIndex;
                value.location = getLocation(groupIndex, voidInstances ? -1 : instanceIndex);
                if (hasLocation2) value.location2 = getLocation2(groupIndex, voidInstances ? -1 : instanceIndex);
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
            value.location2 = NullLocation;
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

export const EmptyLocationIterator: LocationIterator = {
    get hasNext() { return false; },
    get isNextNewInstance() { return false; },
    groupCount: 0,
    instanceCount: 0,
    count: 0,
    stride: 0,
    nonInstanceable: false,
    hasLocation2: false,
    move() {
        return {
            location: NullLocation as Location,
            location2: NullLocation as Location,
            index: 0,
            groupIndex: 0,
            instanceIndex: 0,
            isSecondary: false
        };
    },
    reset() {},
    skipInstance() {},
    voidInstances() {}
};

//

/** A position Location */
export interface PositionLocation {
    readonly kind: 'position-location',
    readonly position: Vec3,
    /** Normal vector at the position (used for surface coloring) */
    readonly normal: Vec3
}
export function PositionLocation(position?: Vec3, normal?: Vec3): PositionLocation {
    return {
        kind: 'position-location',
        position: position ? Vec3.clone(position) : Vec3(),
        normal: normal ? Vec3.clone(normal) : Vec3()
    };
}
export function isPositionLocation(x: any): x is PositionLocation {
    return !!x && x.kind === 'position-location';
}