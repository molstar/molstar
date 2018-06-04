/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { TextureImage, createTextureImage } from 'mol-gl/renderable/util';
import { Color } from 'mol-util/color';
import VertexMap from '../shape/vertex-map';
import { Vec2, Vec3 } from 'mol-math/linear-algebra';

export type ColorType = 'uniform' | 'attribute' | 'instance' | 'element' | 'elementInstance'

export type ColorData = {
    uColor: ValueCell<Vec3>,
    aColor: ValueCell<Float32Array>,
    tColor: ValueCell<TextureImage>,
    uColorTexSize: ValueCell<Vec2>,
    dColorType: ValueCell<string>,
}

const emptyColorTexture = { array: new Uint8Array(3), width: 1, height: 1 }
function createEmptyColorTexture() {
    return {
        tColor: ValueCell.create(emptyColorTexture),
        uColorTexSize: ValueCell.create(Vec2.create(1, 1))
    }
}

export interface UniformColorProps {
    value: Color
}

/** Creates color uniform */
export function createUniformColor(props: UniformColorProps, colorData?: ColorData): ColorData {
    if (colorData) {
        ValueCell.update(colorData.uColor, Color.toRgbNormalized(props.value) as Vec3)
        if (colorData.dColorType.ref.value !== 'uniform') {
            ValueCell.update(colorData.dColorType, 'uniform')
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Color.toRgbNormalized(props.value) as Vec3),
            aColor: ValueCell.create(new Float32Array(0)),
            ...createEmptyColorTexture(),
            dColorType: ValueCell.create('uniform'),
        }
    }
}

export interface AttributeColorProps {
    colorFn: (elementIdx: number) => Color
    vertexMap: VertexMap
}

/** Creates color attribute with color for each element (i.e. shared across instances/units) */
export function createAttributeColor(props: AttributeColorProps, colorData?: ColorData): ColorData {
    const { colorFn, vertexMap } = props
    const { idCount, offsetCount, offsets } = vertexMap
    const colors = new Float32Array(idCount * 3);
    for (let i = 0, il = offsetCount - 1; i < il; ++i) {
        const start = offsets[i]
        const end = offsets[i + 1]
        const hexColor = colorFn(i)
        for (let i = start, il = end; i < il; ++i) {
            Color.toArrayNormalized(hexColor, colors, i * 3)
        }
    }
    if (colorData) {
        ValueCell.update(colorData.aColor, colors)
        if (colorData.dColorType.ref.value !== 'attribute') {
            ValueCell.update(colorData.dColorType, 'attribute')
        }
        return colorData
    } else {
        return {
            uColor: ValueCell.create(Vec3.zero()),
            aColor: ValueCell.create(colors),
            ...createEmptyColorTexture(),
            dColorType: ValueCell.create('attribute'),
        }
    }
}

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

export interface InstanceColorProps {
    colorFn: (instanceIdx: number) => Color
    instanceCount: number
}

/** Creates color texture with color for each instance/unit */
export function createInstanceColor(props: InstanceColorProps, colorData?: ColorData): ColorData {
    const { colorFn, instanceCount} = props
    const colors = colorData && colorData.tColor.ref.value.array.length >= instanceCount * 3 ? colorData.tColor.ref.value : createTextureImage(instanceCount, 3)
    for (let i = 0; i < instanceCount; i++) {
        Color.toArray(colorFn(i), colors.array, i * 3)
    }
    return createTextureColor(colors, 'instance', colorData)
}

export interface ElementColorProps {
    colorFn: (elementIdx: number) => Color
    vertexMap: VertexMap
}

/** Creates color texture with color for each element (i.e. shared across instances/units) */
export function createElementColor(props: ElementColorProps, colorData?: ColorData): ColorData {
    const { colorFn, vertexMap } = props
    const elementCount = vertexMap.offsetCount - 1
    const colors = colorData && colorData.tColor.ref.value.array.length >= elementCount * 3 ? colorData.tColor.ref.value : createTextureImage(elementCount, 3)
    for (let i = 0, il = elementCount; i < il; ++i) {
        Color.toArray(colorFn(i), colors.array, i * 3)
    }
    return createTextureColor(colors, 'element', colorData)
}

export interface ElementInstanceColorProps {
    colorFn: (instanceIdx: number, elementIdx: number) => Color
    instanceCount: number,
    vertexMap: VertexMap
}

/** Creates color texture with color for each element instance (i.e. for each unit) */
export function createElementInstanceColor(props: ElementInstanceColorProps, colorData?: ColorData): ColorData {
    const { colorFn, instanceCount, vertexMap } = props
    const elementCount = vertexMap.offsetCount - 1
    const count = instanceCount * elementCount
    const colors = colorData && colorData.tColor.ref.value.array.length >= count * 3 ? colorData.tColor.ref.value : createTextureImage(count, 3)
    let colorOffset = 0
    for (let i = 0; i < instanceCount; i++) {
        for (let j = 0, jl = elementCount; j < jl; ++j) {
            Color.toArray(colorFn(i, j), colors.array, colorOffset)
            colorOffset += 3
        }
    }
    return createTextureColor(colors, 'elementInstance', colorData)
}

/** Create color attribute or texture, depending on the vertexMap */
export function createAttributeOrElementColor(vertexMap: VertexMap, props: AttributeColorProps, colorData?: ColorData) {
    return vertexMap.idCount < 4 * vertexMap.offsetCount ? createAttributeColor(props, colorData) : createElementColor(props, colorData)
}