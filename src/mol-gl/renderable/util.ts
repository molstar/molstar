/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'
import { ColorData } from 'mol-geo/util/color-data';
import { SizeData } from 'mol-geo/util/size-data';

import { Attributes, AttributesData, AttributesBuffers } from '../renderable'
import Attribute from '../attribute'
import { ShaderDefines } from '../shaders';

export type ReglUniforms = { [k: string]: REGL.Uniform | REGL.Texture }
export type ReglAttributes = { [k: string]: REGL.AttributeConfig }

export function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n * itemSize)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n * itemSize / width) : 0
    return { width, height, length: width * height * itemSize }
}

export interface Texture extends Uint8Array {
    width: number,
    height: number
}

export function createColorTexture (n: number): Texture {
    const colorTexInfo = calculateTextureInfo(n, 3)
    const colorTexture = new Uint8Array(colorTexInfo.length)
    return Object.assign(colorTexture, {
        width: colorTexInfo.width,
        height: colorTexInfo.height
    })
}

export function createTransformAttributes (regl: REGL.Regl, transform: ValueCell<Float32Array>, count: number) {
    const size = 4
    const divisor = 1
    const bpe = transform.ref.value.BYTES_PER_ELEMENT
    const stride = 16 * bpe
    return {
        transformColumn0: Attribute.create(regl, transform, count, { size, divisor, offset: 0, stride }),
        transformColumn1: Attribute.create(regl, transform, count, { size, divisor, offset: 4 * bpe, stride }),
        transformColumn2: Attribute.create(regl, transform, count, { size, divisor, offset: 8 * bpe, stride }),
        transformColumn3: Attribute.create(regl, transform, count, { size, divisor, offset: 12 * bpe, stride })
    }
}

export function createColorUniforms (regl: REGL.Regl, color: ValueCell<Texture>) {
    const colorTex = regl.texture({
        width: color.ref.value.width,
        height: color.ref.value.height,
        format: 'rgb',
        type: 'uint8',
        wrapS: 'clamp',
        wrapT: 'clamp',
        data: color.ref.value
    })
    return {
        colorTex,
        colorTexSize: [ color.ref.value.width, color.ref.value.height ]
    }
}

export function getColorDefines(color: ColorData) {
    const defines: ShaderDefines = {}
    switch (color.type) {
        case 'uniform': defines.UNIFORM_COLOR = ''; break;
        case 'attribute': defines.ATTRIBUTE_COLOR = ''; break;
        case 'element': defines.ELEMENT_COLOR = ''; break;
        case 'instance': defines.INSTANCE_COLOR = ''; break;
        case 'element-instance': defines.ELEMENT_INSTANCE_COLOR = ''; break;
    }
    return defines
}

export function getSizeDefines(size: SizeData) {
    const defines: ShaderDefines = {}
    switch (size.type) {
        case 'uniform': defines.UNIFORM_SIZE = ''; break;
        case 'attribute': defines.ATTRIBUTE_SIZE = ''; break;
    }
    return defines
}

export function getBuffers<T extends AttributesData>(attributes: Attributes<T>): AttributesBuffers<T> {
    const buffers: AttributesBuffers<any> = {}
    for (const k of Object.keys(attributes)) {
        buffers[k] = attributes[k].buffer
    }
    return buffers as AttributesBuffers<T>
}

export function fillSerial<T extends Helpers.NumberArray> (array: T) {
    const n = array.length
    for (let i = 0; i < n; ++i) array[ i ] = i
    return array
}

interface BaseProps {
    objectId: number,
    instanceCount: number,
    elementCount: number,
    positionCount: number,

    position: ValueCell<Float32Array>
    normal?: ValueCell<Float32Array>
    id: ValueCell<Float32Array>
    transform: ValueCell<Float32Array>

    size?: SizeData
    color: ColorData
}

export function createBaseUniforms(regl: REGL.Regl, props: BaseProps): ReglUniforms {
    const { objectId, instanceCount, elementCount, color, size } = props
    const uniforms = { objectId, instanceCount, elementCount }
    if (color.type === 'instance' || color.type === 'element' || color.type === 'element-instance') {
        Object.assign(uniforms, createColorUniforms(regl, color.value))
    } else if (color.type === 'uniform') {
        Object.assign(uniforms, { color: color.value })
    }
    if (size && size.type === 'uniform') {
        Object.assign(uniforms, { size: size.value })
    }
    return uniforms
}

export function createBaseAttributes(regl: REGL.Regl, props: BaseProps): ReglAttributes {
    const { instanceCount, positionCount, position, color, id, normal, size, transform } = props
    const instanceId = ValueCell.create(fillSerial(new Float32Array(instanceCount)))
    const attributes = getBuffers({
        instanceId: Attribute.create(regl, instanceId, instanceCount, { size: 1, divisor: 1 }),
        position: Attribute.create(regl, position, positionCount, { size: 3 }),
        elementId: Attribute.create(regl, id, positionCount, { size: 1 }),
        ...createTransformAttributes(regl, transform, instanceCount)
    })
    if (normal) {
        attributes.normal = Attribute.create(regl, normal as any, positionCount, { size: 3 }).buffer
    }
    if (color.type === 'attribute') {
        attributes.color = Attribute.create(regl, color.value, positionCount, { size: 3 }).buffer
    }
    if (size && size.type === 'attribute') {
        attributes.size = Attribute.create(regl, size.value, positionCount, { size: 1 }).buffer
    }
    return attributes
}

export function createBaseDefines(regl: REGL.Regl, props: BaseProps): ShaderDefines {
    return {
        ...getColorDefines(props.color),
        ...(props.size ? getSizeDefines(props.size) : undefined)
    }
}

export function destroyAttributes(attributes: ReglAttributes) {
    for (const k in attributes) {
        const buffer = attributes[k].buffer
        if (buffer) {
            buffer.destroy()
        }
    }
}

export function destroyUniforms(uniforms: ReglUniforms) {
    for (const k in uniforms) {
        const uniform = uniforms[k]
        if ((uniform as any).destroy) {
            (uniform as any).destroy()
        }
    }
}