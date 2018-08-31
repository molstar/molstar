/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';
import { Color } from 'mol-util/color';
import { Vec2, Vec3 } from 'mol-math/linear-algebra';
import { LocationIterator } from './location-iterator';
import { NullLocation } from 'mol-model/location';
import { LocationColor, ColorType } from 'mol-view/theme/color';

export type ColorData = {
    uColor: ValueCell<Vec3>,
    aColor: ValueCell<Float32Array>,
    tColor: ValueCell<TextureImage>,
    uColorTexDim: ValueCell<Vec2>,
    dColorType: ValueCell<string>,
}

const emptyColorTexture = { array: new Uint8Array(3), width: 1, height: 1 }
function createEmptyColorTexture() {
    return {
        tColor: ValueCell.create(emptyColorTexture),
        uColorTexDim: ValueCell.create(Vec2.create(1, 1))
    }
}

export function createValueColor(value: Color, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.uColor, Color.toRgbNormalized(value) as Vec3)
        if (colorData.dColorType.ref.value !== 'uniform') {
            ValueCell.update(colorData.dColorType, 'uniform')
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Color.toRgbNormalized(value) as Vec3),
            aColor: ValueCell.create(new Float32Array(0)),
            ...createEmptyColorTexture(),
            dColorType: ValueCell.create('uniform'),
        }
    }
}

/** Creates color uniform */
export function createUniformColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    return createValueColor(color(NullLocation, false), colorData)
}

export function createTextureColor(colors: TextureImage, type: ColorType, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.tColor, colors)
        ValueCell.update(colorData.uColorTexDim, Vec2.create(colors.width, colors.height))
        if (colorData.dColorType.ref.value !== type) {
            ValueCell.update(colorData.dColorType, type)
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Vec3.zero()),
            aColor: ValueCell.create(new Float32Array(0)),
            tColor: ValueCell.create(colors),
            uColorTexDim: ValueCell.create(Vec2.create(colors.width, colors.height)),
            dColorType: ValueCell.create(type),
        }
    }
}

/** Creates color texture with color for each instance/unit */
export function createInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { instanceCount} = locationIt
    const colors = colorData && colorData.tColor.ref.value.array.length >= instanceCount * 3 ? colorData.tColor.ref.value : createTextureImage(instanceCount, 3)
    while (locationIt.hasNext) {
        const { location, isSecondary, instanceIndex } = locationIt.move()
        Color.toArray(color(location, isSecondary), colors.array, instanceIndex * 3)
        locationIt.skipInstance()
    }
    return createTextureColor(colors, 'instance', colorData)
}

/** Creates color texture with color for each group (i.e. shared across instances/units) */
export function createGroupColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount } = locationIt
    const colors = colorData && colorData.tColor.ref.value.array.length >= groupCount * 3 ? colorData.tColor.ref.value : createTextureImage(groupCount, 3)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const { location, isSecondary, groupIndex } = locationIt.move()
        Color.toArray(color(location, isSecondary), colors.array, groupIndex * 3)
    }
    return createTextureColor(colors, 'group', colorData)
}

/** Creates color texture with color for each group in each instance (i.e. for each unit) */
export function createGroupInstanceColor(locationIt: LocationIterator, color: LocationColor, colorData?: ColorData): ColorData {
    const { groupCount, instanceCount } = locationIt
    const count = instanceCount * groupCount
    console.log(count, instanceCount, groupCount)
    const colors = colorData && colorData.tColor.ref.value.array.length >= count * 3 ? colorData.tColor.ref.value : createTextureImage(count, 3)
    console.log(colors.array.length / 3, count)
    while (locationIt.hasNext) {
        const { location, isSecondary, index } = locationIt.move()
        Color.toArray(color(location, isSecondary), colors.array, index * 3)
    }
    return createTextureColor(colors, 'groupInstance', colorData)
}