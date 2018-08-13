/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';
import { Color } from 'mol-util/color';
import { Vec2, Vec3 } from 'mol-math/linear-algebra';
import { LocationIterator } from '../representation/structure/visual/util/location-iterator';
import { Location, NullLocation } from 'mol-model/location';

export type ColorType = 'uniform' | 'attribute' | 'instance' | 'element' | 'elementInstance'

export type ColorData = {
    uColor: ValueCell<Vec3>,
    aColor: ValueCell<Float32Array>,
    tColor: ValueCell<TextureImage>,
    uColorTexSize: ValueCell<Vec2>,
    dColorType: ValueCell<string>,
}

export type LocationColor = (location: Location, isSecondary: boolean) => Color

const emptyColorTexture = { array: new Uint8Array(3), width: 1, height: 1 }
function createEmptyColorTexture() {
    return {
        tColor: ValueCell.create(emptyColorTexture),
        uColorTexSize: ValueCell.create(Vec2.create(1, 1))
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
export function createUniformColor(locationIt: LocationIterator, colorFn: LocationColor, colorData?: ColorData): ColorData {
    return createValueColor(colorFn(NullLocation, false), colorData)
}

// export interface AttributeColorProps {
//     colorFn: (elementIdx: number) => Color
//     vertexMap: VertexMap
// }

// /** Creates color attribute with color for each element (i.e. shared across instances/units) */
// export function createAttributeColor(props: AttributeColorProps, colorData?: ColorData): ColorData {
//     const { colorFn, vertexMap } = props
//     const { idCount, offsetCount, offsets } = vertexMap
//     const colors = new Float32Array(idCount * 3);
//     for (let i = 0, il = offsetCount - 1; i < il; ++i) {
//         const start = offsets[i]
//         const end = offsets[i + 1]
//         const hexColor = colorFn(i)
//         for (let i = start, il = end; i < il; ++i) {
//             Color.toArrayNormalized(hexColor, colors, i * 3)
//         }
//     }
//     if (colorData) {
//         ValueCell.update(colorData.aColor, colors)
//         if (colorData.dColorType.ref.value !== 'attribute') {
//             ValueCell.update(colorData.dColorType, 'attribute')
//         }
//         return colorData
//     } else {
//         return {
//             uColor: ValueCell.create(Vec3.zero()),
//             aColor: ValueCell.create(colors),
//             ...createEmptyColorTexture(),
//             dColorType: ValueCell.create('attribute'),
//         }
//     }
// }

export function createTextureColor(colors: TextureImage, type: ColorType, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.tColor, colors)
        ValueCell.update(colorData.uColorTexSize, Vec2.create(colors.width, colors.height))
        if (colorData.dColorType.ref.value !== type) {
            ValueCell.update(colorData.dColorType, type)
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Vec3.zero()),
            aColor: ValueCell.create(new Float32Array(0)),
            tColor: ValueCell.create(colors),
            uColorTexSize: ValueCell.create(Vec2.create(colors.width, colors.height)),
            dColorType: ValueCell.create(type),
        }
    }
}

/** Creates color texture with color for each instance/unit */
export function createInstanceColor(locationIt: LocationIterator, colorFn: LocationColor, colorData?: ColorData): ColorData {
    const { instanceCount} = locationIt
    const colors = colorData && colorData.tColor.ref.value.array.length >= instanceCount * 3 ? colorData.tColor.ref.value : createTextureImage(instanceCount, 3)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move()
        Color.toArray(colorFn(v.location, v.isSecondary), colors.array, v.index * 3)
        locationIt.skipInstance()
    }
    return createTextureColor(colors, 'instance', colorData)
}

/** Creates color texture with color for each element (i.e. shared across instances/units) */
export function createElementColor(locationIt: LocationIterator, colorFn: LocationColor, colorData?: ColorData): ColorData {
    const { elementCount } = locationIt
    const colors = colorData && colorData.tColor.ref.value.array.length >= elementCount * 3 ? colorData.tColor.ref.value : createTextureImage(elementCount, 3)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move()
        // console.log(v)
        Color.toArray(colorFn(v.location, v.isSecondary), colors.array, v.elementIndex * 3)
    }
    return createTextureColor(colors, 'element', colorData)
}

/** Creates color texture with color for each element instance (i.e. for each unit) */
export function createElementInstanceColor(locationIt: LocationIterator, colorFn: LocationColor, colorData?: ColorData): ColorData {
    const { elementCount, instanceCount } = locationIt
    const count = instanceCount * elementCount
    const colors = colorData && colorData.tColor.ref.value.array.length >= count * 3 ? colorData.tColor.ref.value : createTextureImage(count, 3)
    while (locationIt.hasNext && !locationIt.isNextNewInstance) {
        const v = locationIt.move()
        Color.toArray(colorFn(v.location, v.isSecondary), colors.array, v.index * 3)
    }
    return createTextureColor(colors, 'elementInstance', colorData)
}