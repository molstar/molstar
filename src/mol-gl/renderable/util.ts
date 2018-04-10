/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'

import { ColorTexture } from '../util';
import { Attributes, AttributesData, AttributesBuffers } from '../renderable'
import Attribute from '../attribute'

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

export function createColorUniforms (regl: REGL.Regl, color: ValueCell<ColorTexture>) {
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