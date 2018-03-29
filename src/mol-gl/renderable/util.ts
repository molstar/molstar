/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');

import Attribute from '../attribute'
import { Attributes, AttributesData, AttributesBuffers } from '../renderable'

export function getData<T extends AttributesData>(attributes: Attributes<T>): T {
    const data: AttributesData = {}
    for (const k of Object.keys(attributes)) {
        data[k] = attributes[k].data
    }
    return data as T
}

export function getBuffers<T extends AttributesData>(attributes: Attributes<T>): AttributesBuffers<T> {
    const buffers: AttributesBuffers<any> = {}
    for (const k of Object.keys(attributes)) {
        buffers[k] = attributes[k].buffer
    }
    return buffers as AttributesBuffers<T>
}

export function createAttributes<T extends AttributesData>(regl: REGL.Regl, data: T): Attributes<T> {
    const attributes: Attributes<any> = {}
    for (const k of Object.keys(data)) {
        attributes[k] = Attribute.create(regl, data[k])
    }
    return attributes as Attributes<T>
}