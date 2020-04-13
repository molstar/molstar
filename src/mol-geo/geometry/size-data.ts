/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { Vec2 } from '../../mol-math/linear-algebra';
import { TextureImage, createTextureImage } from '../../mol-gl/renderable/util';
import { LocationIterator } from '../util/location-iterator';
import { Location, NullLocation } from '../../mol-model/location';
import { SizeTheme } from '../../mol-theme/size';
import { Geometry } from './geometry';
import { encodeFloatLog, decodeFloatLog } from '../../mol-util/float-packing';

export type SizeType = 'uniform' | 'instance' | 'group' | 'groupInstance'

export type SizeData = {
    uSize: ValueCell<number>,
    tSize: ValueCell<TextureImage<Uint8Array>>,
    uSizeTexDim: ValueCell<Vec2>,
    dSizeType: ValueCell<string>,
}

export function createSizes(locationIt: LocationIterator, sizeTheme: SizeTheme<any>, sizeData?: SizeData): SizeData {
    switch (Geometry.getGranularity(locationIt, sizeTheme.granularity)) {
        case 'uniform': return createUniformSize(locationIt, sizeTheme.size, sizeData);
        case 'group': return createGroupSize(locationIt, sizeTheme.size, sizeData);
        case 'groupInstance': return createGroupInstanceSize(locationIt, sizeTheme.size, sizeData);
        case 'instance': return createInstanceSize(locationIt, sizeTheme.size, sizeData);
    }
}

export function getMaxSize(sizeData: SizeData): number {
    const type = sizeData.dSizeType.ref.value as SizeType;
    switch (type) {
        case 'uniform':
            return sizeData.uSize.ref.value;
        case 'instance':
        case 'group':
        case 'groupInstance':
            let maxSize = 0;
            const array = sizeData.tSize.ref.value.array;
            for (let i = 0, il = array.length; i < il; ++i) {
                const value = decodeFloatLog(array[i] / 255);
                if (maxSize < value) maxSize = value;
            }
            return maxSize;
    }
}

export type LocationSize = (location: Location) => number

const emptySizeTexture = { array: new Uint8Array(1), width: 1, height: 1 };
function createEmptySizeTexture() {
    return {
        tSize: ValueCell.create(emptySizeTexture),
        uSizeTexDim: ValueCell.create(Vec2.create(1, 1))
    };
}

export function createValueSize(value: number, sizeData?: SizeData): SizeData {
    if (sizeData) {
        ValueCell.update(sizeData.uSize, value);
        if (sizeData.dSizeType.ref.value !== 'uniform') {
            ValueCell.update(sizeData.dSizeType, 'uniform');
        }
        return sizeData;
    } else {
        return {
            uSize: ValueCell.create(value),
            ...createEmptySizeTexture(),
            dSizeType: ValueCell.create('uniform'),
        };
    }
}

/** Creates size uniform */
export function createUniformSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    return createValueSize(sizeFn(NullLocation), sizeData);
}

export function createTextureSize(sizes: TextureImage<Uint8Array>, type: SizeType, sizeData?: SizeData): SizeData {
    if (sizeData) {
        ValueCell.update(sizeData.tSize, sizes);
        ValueCell.update(sizeData.uSizeTexDim, Vec2.create(sizes.width, sizes.height));
        if (sizeData.dSizeType.ref.value !== type) {
            ValueCell.update(sizeData.dSizeType, type);
        }
        return sizeData;
    } else {
        return {
            uSize: ValueCell.create(0),
            tSize: ValueCell.create(sizes),
            uSizeTexDim: ValueCell.create(Vec2.create(sizes.width, sizes.height)),
            dSizeType: ValueCell.create(type),
        };
    }
}

/** Creates size texture with size for each instance/unit */
export function createInstanceSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    const { instanceCount} = locationIt;
    const sizes = createTextureImage(Math.max(1, instanceCount), 1, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move();
        sizes.array[v.instanceIndex] = encodeFloatLog(sizeFn(v.location)) * 255;
        locationIt.skipInstance();
    }
    return createTextureSize(sizes, 'instance', sizeData);
}

/** Creates size texture with size for each group (i.e. shared across instances/units) */
export function createGroupSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    const { groupCount } = locationIt;
    const sizes = createTextureImage(Math.max(1, groupCount), 1, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move();
        sizes.array[v.groupIndex] = encodeFloatLog(sizeFn(v.location)) * 255;
    }
    return createTextureSize(sizes, 'group', sizeData);
}

/** Creates size texture with size for each group in each instance (i.e. for each unit) */
export function createGroupInstanceSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    const { groupCount, instanceCount } = locationIt;
    const count = instanceCount * groupCount;
    const sizes = createTextureImage(Math.max(1, count), 1, Uint8Array, sizeData && sizeData.tSize.ref.value.array);
    locationIt.reset();
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move();
        sizes.array[v.index] = encodeFloatLog(sizeFn(v.location)) * 255;
    }
    return createTextureSize(sizes, 'groupInstance', sizeData);
}