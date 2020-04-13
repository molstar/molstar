/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Iterator } from '../../mol-data';
import { NullLocation, Location } from '../../mol-model/location';

export interface LocationValue {
    location: Location
    index: number
    groupIndex: number
    instanceIndex: number
    isSecondary: boolean
}

export const NullLocationValue: LocationValue = {
    location: NullLocation,
    index: 0,
    groupIndex: 0,
    instanceIndex: 0,
    isSecondary: false
};

export interface LocationIterator extends Iterator<LocationValue> {
    readonly hasNext: boolean
    readonly isNextNewInstance: boolean
    readonly groupCount: number
    readonly instanceCount: number
    readonly count: number
    /** If true, may have multiple units per instance; if false one unit per instance */
    readonly isComplex: boolean
    move(): LocationValue
    reset(): void
    skipInstance(): void
}

type LocationGetter = (groupIndex: number, instanceIndex: number) => Location
type IsSecondaryGetter = (groupIndex: number, instanceIndex: number) => boolean

export function LocationIterator(groupCount: number, instanceCount: number, getLocation: LocationGetter, isComplex = false, isSecondary: IsSecondaryGetter = () => false): LocationIterator {
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

    return {
        get hasNext () { return hasNext; },
        get isNextNewInstance () { return isNextNewInstance; },
        groupCount,
        instanceCount,
        count: groupCount * instanceCount,
        isComplex,
        move() {
            if (hasNext) {
                value.groupIndex = groupIndex;
                value.instanceIndex = instanceIndex;
                value.index = instanceIndex * groupCount + groupIndex;
                value.location = getLocation(groupIndex, instanceIndex);
                value.isSecondary = isSecondary(groupIndex, instanceIndex);
                ++groupIndex;
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
        },
        skipInstance() {
            if (hasNext && value.instanceIndex === instanceIndex) {
                ++instanceIndex;
                groupIndex = 0;
                hasNext = instanceIndex < instanceCount;
            }
        }
    };
}