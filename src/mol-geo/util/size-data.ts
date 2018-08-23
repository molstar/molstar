/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { Vec2 } from 'mol-math/linear-algebra';
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';
import { LocationIterator } from './location-iterator';
import { Location, NullLocation } from 'mol-model/location';

export type SizeType = 'uniform' | 'instance' | 'group' | 'groupInstance'

export type SizeData = {
    uSize: ValueCell<number>,
    aSize: ValueCell<Float32Array>,
    tSize: ValueCell<TextureImage>,
    uSizeTexDim: ValueCell<Vec2>,
    dSizeType: ValueCell<string>,
}

export type LocationSize = (location: Location) => number

const emptySizeTexture = { array: new Uint8Array(1), width: 1, height: 1 }
function createEmptySizeTexture() {
    return {
        tSize: ValueCell.create(emptySizeTexture),
        uSizeTexDim: ValueCell.create(Vec2.create(1, 1))
    }
}

export function createValueSize(value: number, sizeData?: SizeData): SizeData {
    if (sizeData) {
        ValueCell.update(sizeData.uSize, value)
        if (sizeData.dSizeType.ref.value !== 'uniform') {
            ValueCell.update(sizeData.dSizeType, 'uniform')
        }
        return sizeData
    } else {
        return {
            uSize: ValueCell.create(value),
            aSize: ValueCell.create(new Float32Array(0)),
            ...createEmptySizeTexture(),
            dSizeType: ValueCell.create('uniform'),
        }
    }
}

/** Creates size uniform */
export function createUniformSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    return createValueSize(sizeFn(NullLocation), sizeData)
}

export function createTextureSize(sizes: TextureImage, type: SizeType, sizeData?: SizeData): SizeData {
    if (sizeData) {
        ValueCell.update(sizeData.tSize, sizes)
        ValueCell.update(sizeData.uSizeTexDim, Vec2.create(sizes.width, sizes.height))
        if (sizeData.dSizeType.ref.value !== type) {
            ValueCell.update(sizeData.dSizeType, type)
        }
        return sizeData
    } else {
        return {
            uSize: ValueCell.create(0),
            aSize: ValueCell.create(new Float32Array(0)),
            tSize: ValueCell.create(sizes),
            uSizeTexDim: ValueCell.create(Vec2.create(sizes.width, sizes.height)),
            dSizeType: ValueCell.create(type),
        }
    }
}

/** Creates size texture with size for each instance/unit */
export function createInstanceSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    const { instanceCount} = locationIt
    const sizes = sizeData && sizeData.tSize.ref.value.array.length >= instanceCount ? sizeData.tSize.ref.value : createTextureImage(instanceCount, 1)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move()
        sizes.array[v.instanceIndex] = sizeFn(v.location)
        locationIt.skipInstance()
    }
    return createTextureSize(sizes, 'instance', sizeData)
}

/** Creates size texture with size for each group (i.e. shared across instances/units) */
export function createGroupSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    const { groupCount } = locationIt
    const sizes = sizeData && sizeData.tSize.ref.value.array.length >= groupCount ? sizeData.tSize.ref.value : createTextureImage(groupCount, 1)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move()
        sizes.array[v.groupIndex] = sizeFn(v.location)
    }
    return createTextureSize(sizes, 'group', sizeData)
}

/** Creates size texture with size for each group in each instance (i.e. for each unit) */
export function createGroupInstanceSize(locationIt: LocationIterator, sizeFn: LocationSize, sizeData?: SizeData): SizeData {
    const { groupCount, instanceCount } = locationIt
    const count = instanceCount * groupCount
    const sizes = sizeData && sizeData.tSize.ref.value.array.length >= count ? sizeData.tSize.ref.value : createTextureImage(count, 1)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move()
        sizes.array[v.index] = sizeFn(v.location)
    }
    return createTextureSize(sizes, 'groupInstance', sizeData)
}