/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { Renderable } from '../renderable'
import { getBuffers } from './util'
import Attribute from '../attribute';

import { MeshShaders } from '../shaders'

type Mesh = 'mesh'

type Uniforms = { [k: string]: REGL.Uniform | REGL.Texture }

export function fillSerial<T extends Helpers.NumberArray> (array: T) {
    const n = array.length
    for (let i = 0; i < n; ++i) array[ i ] = i
    return array
}

namespace Mesh {
    export type DataType = {
        position: { type: Float32Array, itemSize: 3 }
        normal: { type: Float32Array, itemSize: 3 }
        transformColumn0: { type: Float32Array, itemSize: 4 }
        transformColumn1: { type: Float32Array, itemSize: 4 }
        transformColumn2: { type: Float32Array, itemSize: 4 }
        transformColumn3: { type: Float32Array, itemSize: 4 }
    }
    export type Data = { [K in keyof DataType]: DataType[K]['type'] }
    export type Attributes = { [K in keyof Data]: Attribute<Data[K]> }

    export function create(regl: REGL.Regl, attributes: Attributes, uniforms: Uniforms, elements?: Helpers.UintArray): Renderable<Data> {
        // console.log('mesh', {
        //     count: attributes.position.getCount(),
        //     instances: attributes.transformColumn0.getCount(),
        //     attributes,
        //     uniforms
        // })
        const instanceCount = attributes.transformColumn0.getCount()
        const instanceId = fillSerial(new Float32Array(instanceCount))
        // console.log(instanceId)
        const command = regl({
            ...MeshShaders,
            uniforms: {
                objectId: uniforms.objectId || 0,
                instanceCount,
                ...uniforms
            },
            attributes: getBuffers({
                instanceId: Attribute.create(regl, instanceId, { size: 1, divisor: 1 }),
                ...attributes
            }),
            elements: elements && regl.elements({
                data: new Uint16Array(elements),
                primitive: 'triangles',
                // type: 'uint16',
                // count: elements.length / 3,
                // length: elements.length * 2
            }),
            count: elements ? elements.length : attributes.position.getCount(),
            instances: instanceCount,
            primitive: 'triangles'
        })
        return {
            draw: () => command(),
        }
    }
}

export default Mesh