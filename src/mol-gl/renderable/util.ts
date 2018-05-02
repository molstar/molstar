/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { ColorData } from 'mol-geo/util/color-data';
import { SizeData } from 'mol-geo/util/size-data';

import { ShaderDefines } from '../shader-code';
import { UniformDefs, UniformValues } from '../webgl/uniform';
import { AttributeDefs, AttributeValues } from '../webgl/buffer';
import { Vec3, Vec2 } from 'mol-math/linear-algebra';
import { TextureDefs, TextureValues } from '../webgl/texture';

export function calculateTextureInfo (n: number, itemSize: number) {
    const sqN = Math.sqrt(n * itemSize)
    let width = Math.ceil(sqN)
    width = width + (itemSize - (width % itemSize)) % itemSize
    const height = width > 0 ? Math.ceil(n * itemSize / width) : 0
    return { width, height, length: width * height * itemSize }
}

export interface TextureImage {
    array: Uint8Array
    width: number
    height: number
}

export function createColorTexture (n: number): TextureImage {
    const { length, width, height } = calculateTextureInfo(n, 3)
    return { array: new Uint8Array(length), width, height }
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
    normal?: ValueCell<Float32Array | undefined>
    id: ValueCell<Float32Array>
    transform: ValueCell<Float32Array>

    size?: SizeData
    color: ColorData
}

export function getBaseUniformDefs(props: BaseProps) {
    const uniformDefs: UniformDefs = {
        model: 'm4',
        view: 'm4',
        projection: 'm4',

        pixelRatio: 'f',
        viewportHeight: 'f',

        // light_position: 'v3',
        light_color: 'v3',
        light_ambient: 'v3',

        objectId: 'i',
        instanceCount: 'i',
        elementCount: 'i'
    }
    const color = props.color
    if (color.type === 'instance' || color.type === 'element' || color.type === 'element-instance') {
        // uniformDefs.colorTex = 't2'
        uniformDefs.colorTexSize = 'v2'
    } else if (color.type === 'uniform') {
        uniformDefs.color = 'v3'
    }
    const size = props.size
    if (size && size.type === 'uniform') {
        uniformDefs.size = 'f'
    }
    return uniformDefs
}

export function getBaseUniformValues(props: BaseProps) {
    const { objectId, instanceCount, elementCount } = props
    const uniformValues: UniformValues = {
        objectId, instanceCount, elementCount
    }
    const color = props.color
    if (color.type === 'instance' || color.type === 'element' || color.type === 'element-instance') {
        const { width, height } = color.data.ref.value
        // uniformValues.colorTex = color.value.ref.value.array
        uniformValues.colorTexSize = Vec2.create(width, height)
    } else if (color.type === 'uniform') {
        uniformValues.color = color.data as Vec3
    }
    const size = props.size
    if (size && size.type === 'uniform') {
        uniformValues.size = size.value
    }
    return uniformValues
}

export function getBaseAttributeDefs(props: BaseProps) {
    const attributeDefs: AttributeDefs = {
        instanceId: { kind: 'float32', itemSize: 1, divisor: 1 },
        position: { kind: 'float32', itemSize: 3, divisor: 0 },
        elementId: { kind: 'float32', itemSize: 1, divisor: 0 },
        transform: { kind: 'float32', itemSize: 16, divisor: 1 },
    }
    if (props.normal && props.normal.ref.value) {
        attributeDefs.normal = { kind: 'float32', itemSize: 3, divisor: 0 }
    }
    const color = props.color
    if (color.type === 'attribute') {
        attributeDefs.color = { kind: 'float32', itemSize: 3, divisor: 0 }
    }
    const size = props.size
    if (size && size.type === 'attribute') {
        attributeDefs.size = { kind: 'float32', itemSize: 1, divisor: 0 }
    }
    return attributeDefs
}

export function getBaseAttributeValues(props: BaseProps) {
    const { instanceCount, position, id, normal, transform } = props
    const instanceId = ValueCell.create(fillSerial(new Float32Array(instanceCount)))
    const attributeValues: AttributeValues = {
        instanceId: instanceId.ref.value,
        position: position.ref.value,
        elementId: id.ref.value,
        transform: transform.ref.value
    }
    if (normal && normal.ref.value) {
        attributeValues.normal = normal.ref.value
    }
    const color = props.color
    if (color.type === 'attribute') {
        attributeValues.color = color.data.ref.value
    }
    const size = props.size
    if (size && size.type === 'attribute') {
        attributeValues.size = size.value.ref.value
    }
    return attributeValues
}

export function getBaseTextureDefs(props: BaseProps) {
    const textureDefs: TextureDefs = {}
    const color = props.color
    if (color.type === 'instance' || color.type === 'element' || color.type === 'element-instance') {
        textureDefs.colorTex = true
    }
    return textureDefs
}

export function getBaseTextureValues(props: BaseProps) {
    const textureValues: TextureValues = {}
    const color = props.color
    if (color.type === 'instance' || color.type === 'element' || color.type === 'element-instance') {
        textureValues.colorTex = color.data.ref.value
    }
    return textureValues
}

export function getBaseDefines(props: BaseProps): ShaderDefines {
    return {
        ...getColorDefines(props.color),
        ...(props.size ? getSizeDefines(props.size) : undefined)
    }
}

export function getBaseDefs(props: BaseProps) {
    return {
        uniformDefs: getBaseUniformDefs(props),
        attributeDefs: getBaseAttributeDefs(props),
        textureDefs: getBaseTextureDefs(props),
    }
}

export function getBaseValues(props: BaseProps) {
    return {
        uniformValues: getBaseUniformValues(props),
        attributeValues: getBaseAttributeValues(props),
        textureValues: getBaseTextureValues(props),
    }
}
