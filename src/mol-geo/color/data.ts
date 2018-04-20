/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util';
import Color from './color';
import { Mesh } from '../shape/mesh';

function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n * itemSize)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n * itemSize / width) : 0
    return { width, height, length: width * height * itemSize }
}

interface ColorTexture extends Uint8Array {
    width: number,
    height: number
}

function createColorTexture (n: number): ColorTexture {
    const colorTexInfo = calculateTextureInfo(n, 3)
    const colorTexture = new Uint8Array(colorTexInfo.length)
    return Object.assign(colorTexture, {
        width: colorTexInfo.width,
        height: colorTexInfo.height
    })
}

export type UniformColor = { type: 'uniform', value: number[] }
export type AttributeColor = { type: 'attribute', value: ValueCell<Float32Array> }
export type InstanceColor = { type: 'instance', value: ValueCell<ColorTexture> }
export type ElementColor = { type: 'element', value: ValueCell<ColorTexture> }
export type ElementInstanceColor = { type: 'element-instance', value: ValueCell<ColorTexture> }
export type ColorData = UniformColor | AttributeColor | InstanceColor | ElementColor | ElementInstanceColor

export interface UniformColorProps {
    value: Color
}

/** Creates color uniform */
export function createUniformColor(props: UniformColorProps): UniformColor {
    return { type: 'uniform', value: Color.toRgbNormalized(props.value) }
}

export interface AttributeColorProps {
    colorFn: (elementIdx: number) => Color
    vertexCount: number,
    offsetCount: number,
    offsets: ValueCell<Uint32Array>
}

/** Creates color attribute with color for each element (i.e. shared across indtances/units) */
export function createAttributeColor(props: AttributeColorProps): AttributeColor {
    const { colorFn, vertexCount, offsetCount, offsets} = props
    const colors = new Float32Array(vertexCount * 3);
    const _offsets = offsets.ref.value
    for (let i = 0, il = offsetCount - 1; i < il; ++i) {
        const start = _offsets[i]
        const end = _offsets[i + 1]
        const hexColor = colorFn(i)
        for (let i = start, il = end; i < il; ++i) {
            Color.toArrayNormalized(hexColor, colors, i * 3)
        }
    }
    return { type: 'attribute', value: ValueCell.create(colors) }
}

export interface InstanceColorProps {
    colorFn: (unitIdx: number) => Color
    unitCount: number
}

/** Creates color texture with color for each instance/unit */
export function createInstanceColor(props: InstanceColorProps): InstanceColor {
    const { colorFn, unitCount} = props
    const colors = createColorTexture(unitCount)
    for (let i = 0; i < unitCount; i++) {
        Color.toArray(colorFn(i), colors, i * 3)
    }
    return { type: 'instance', value: ValueCell.create(colors) }
}

export interface ElementColorProps {
    colorFn: (elementIdx: number) => Color
    offsetCount: number,
    offsets: ValueCell<Uint32Array>
}

/** Creates color texture with color for each element (i.e. shared across indtances/units) */
export function createElementColor(props: ElementColorProps): ElementColor {
    const { colorFn, offsetCount } = props
    const elementCount = offsetCount - 1
    const colors = createColorTexture(elementCount)
    for (let i = 0, il = elementCount; i < il; ++i) {
        Color.toArray(colorFn(i), colors, i * 3)
    }
    return { type: 'element', value: ValueCell.create(colors) }
}

export interface ElementInstanceColorProps {
    colorFn: (unitIdx: number, elementIdx: number) => Color
    unitCount: number,
    offsetCount: number,
    offsets: ValueCell<Uint32Array>
}

/** Creates color texture with color for each element instance (i.e. for each unit) */
export function createElementInstanceColor(props: ElementInstanceColorProps): ElementInstanceColor {
    const { colorFn, unitCount, offsetCount } = props
    const elementCount = offsetCount - 1
    const count = unitCount * elementCount
    const colors = createColorTexture(count)
    let colorOffset = 0
    for (let i = 0; i < unitCount; i++) {
        for (let j = 0, jl = elementCount; j < jl; ++j) {
            Color.toArray(colorFn(i, j), colors, colorOffset)
            colorOffset += 3
        }
    }
    return { type: 'element-instance', value: ValueCell.create(colors) }
}

/** Create color attribute or texture, depending on the mesh */
export function createAttributeOrElementColor(mesh: Mesh, props: AttributeColorProps) {
    return mesh.vertexCount < 4 * mesh.offsetCount ? createAttributeColor(props) : createElementColor(props)
}