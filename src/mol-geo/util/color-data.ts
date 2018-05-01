/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import { TextureImage, createColorTexture } from 'mol-gl/renderable/util';
import { Color } from 'mol-util/color';
import VertexMap from '../shape/vertex-map';

export type UniformColor = { type: 'uniform', data: number[] }
export type AttributeColor = { type: 'attribute', data: ValueCell<Float32Array> }
export type InstanceColor = { type: 'instance', data: ValueCell<TextureImage> }
export type ElementColor = { type: 'element', data: ValueCell<TextureImage> }
export type ElementInstanceColor = { type: 'element-instance', data: ValueCell<TextureImage> }
export type ColorData = UniformColor | AttributeColor | InstanceColor | ElementColor | ElementInstanceColor

export interface UniformColorProps {
    value: Color
}

/** Creates color uniform */
export function createUniformColor(props: UniformColorProps): UniformColor {
    return { type: 'uniform', data: Color.toRgbNormalized(props.value) }
}

export interface AttributeColorProps {
    colorFn: (elementIdx: number) => Color
    vertexMap: VertexMap
}

/** Creates color attribute with color for each element (i.e. shared across instances/units) */
export function createAttributeColor(props: AttributeColorProps): AttributeColor {
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
    return { type: 'attribute', data: ValueCell.create(colors) }
}

export interface InstanceColorProps {
    colorFn: (instanceIdx: number) => Color
    instanceCount: number
}

/** Creates color texture with color for each instance/unit */
export function createInstanceColor(props: InstanceColorProps): InstanceColor {
    const { colorFn, instanceCount} = props
    const colors = createColorTexture(instanceCount)
    for (let i = 0; i < instanceCount; i++) {
        Color.toArray(colorFn(i), colors.array, i * 3)
    }
    return { type: 'instance', data: ValueCell.create(colors) }
}

export interface ElementColorProps {
    colorFn: (elementIdx: number) => Color
    vertexMap: VertexMap
}

/** Creates color texture with color for each element (i.e. shared across instances/units) */
export function createElementColor(props: ElementColorProps): ElementColor {
    const { colorFn, vertexMap } = props
    const elementCount = vertexMap.offsetCount - 1
    const colors = createColorTexture(elementCount)
    for (let i = 0, il = elementCount; i < il; ++i) {
        Color.toArray(colorFn(i), colors.array, i * 3)
    }
    return { type: 'element', data: ValueCell.create(colors) }
}

export interface ElementInstanceColorProps {
    colorFn: (instanceIdx: number, elementIdx: number) => Color
    instanceCount: number,
    vertexMap: VertexMap
}

/** Creates color texture with color for each element instance (i.e. for each unit) */
export function createElementInstanceColor(props: ElementInstanceColorProps): ElementInstanceColor {
    const { colorFn, instanceCount, vertexMap } = props
    const elementCount = vertexMap.offsetCount - 1
    const count = instanceCount * elementCount
    const colors = createColorTexture(count)
    let colorOffset = 0
    for (let i = 0; i < instanceCount; i++) {
        for (let j = 0, jl = elementCount; j < jl; ++j) {
            Color.toArray(colorFn(i, j), colors.array, colorOffset)
            colorOffset += 3
        }
    }
    return { type: 'element-instance', data: ValueCell.create(colors) }
}

/** Create color attribute or texture, depending on the vertexMap */
export function createAttributeOrElementColor(vertexMap: VertexMap, props: AttributeColorProps) {
    return vertexMap.idCount < 4 * vertexMap.offsetCount ? createAttributeColor(props) : createElementColor(props)
}