/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { Renderable } from '../renderable'
import { getBuffers } from './util'
import Attribute from '../attribute';
import { PointShaders } from '../shaders'

type Point = 'point'

namespace Point {
    export type DataType = {
        position: { type: Float32Array, itemSize: 3 }
        transformColumn0: { type: Float32Array, itemSize: 4 }
        transformColumn1: { type: Float32Array, itemSize: 4 }
        transformColumn2: { type: Float32Array, itemSize: 4 }
        transformColumn3: { type: Float32Array, itemSize: 4 }
    }
    export type Data = { [K in keyof DataType]: DataType[K]['type'] }
    export type Attributes = { [K in keyof Data]: Attribute<Data[K]> }

    export function create(regl: REGL.Regl, attributes: Attributes): Renderable<Data> {
        console.log('point', {
            count: attributes.position.getCount(),
            instances: attributes.transformColumn0.getCount(),
        }, attributes)
        const command = regl({
            ...PointShaders,
            attributes: getBuffers(attributes),
            count: attributes.position.getCount(),
            instances: attributes.transformColumn0.getCount(),
            primitive: 'points'
        })
        return {
            draw: () => command(),
        }
    }
}

export default Point

// namespace Point {
//     export type DataType = {
//         position: { type: Float32Array, itemSize: 3 }
//     }
//     export type Data = { [K in keyof DataType]: DataType[K]['type'] }
//     export type Attributes = { [K in keyof Data]: Attribute<Data[K]> }

//     export function create(regl: REGL.Regl, dataOrCount: Data | number): Renderable<Data> {
//         let count: number
//         let data: Data
//         if (typeof dataOrCount === 'number') {
//             count = dataOrCount
//             data = {
//                 position: new Float32Array(count * 3)
//             }
//         } else {
//             count = dataOrCount.position.length / 3
//             data = dataOrCount
//         }
//         const attributes = createAttributes(regl, data)
//         const command = regl({
//             vert: pointVert,
//             frag: pointFrag,
//             attributes: getBuffers(attributes),
//             count,
//             primitive: 'points'
//         })
//         return {
//             draw: () => command(),
//             setCount: (newCount: number) => {
//                 for (const k of Object.keys(data)) {
//                     attributes[k as keyof Data].setCount(newCount)
//                 }
//                 count = newCount
//             },
//             getCount: () => count,
//             attributes
//         }
//     }
// }