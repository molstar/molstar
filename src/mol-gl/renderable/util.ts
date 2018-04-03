/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { Attributes, AttributesData, AttributesBuffers } from '../renderable'
import Attribute from '../attribute'

export function createTransformAttributes (regl: REGL.Regl, transform: Float32Array) {
    const size = 4
    const divisor = 1
    const bpe = transform.BYTES_PER_ELEMENT
    const stride = 16 * bpe
    return {
        transformColumn0: Attribute.create(regl, transform, { size, divisor, offset: 0, stride }),
        transformColumn1: Attribute.create(regl, transform, { size, divisor, offset: 4 * bpe, stride }),
        transformColumn2: Attribute.create(regl, transform, { size, divisor, offset: 8 * bpe, stride }),
        transformColumn3: Attribute.create(regl, transform, { size, divisor, offset: 12 * bpe, stride })
    }
}

export function getBuffers<T extends AttributesData>(attributes: Attributes<T>): AttributesBuffers<T> {
    const buffers: AttributesBuffers<any> = {}
    for (const k of Object.keys(attributes)) {
        buffers[k] = attributes[k].buffer
    }
    return buffers as AttributesBuffers<T>
}
