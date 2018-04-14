/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'

import { Renderable } from '../renderable'
import { getBuffers, createTransformAttributes, fillSerial } from './util'
import Attribute from '../attribute';
import { PointShaders } from '../shaders'

type Point = 'point'

namespace Point {
    export type Data = {
        position: ValueCell<Float32Array>
        transform: ValueCell<Float32Array>

        instanceCount: number
        positionCount: number
    }

    export function create(regl: REGL.Regl, data: Data): Renderable {
        const instanceId = ValueCell.create(fillSerial(new Float32Array(data.instanceCount)))
        const command = regl({
            ...PointShaders,
            attributes: getBuffers({
                instanceId: Attribute.create(regl, instanceId, data.instanceCount, { size: 1, divisor: 1 }),
                position: Attribute.create(regl, data.position, data.positionCount, { size: 3 }),
                ...createTransformAttributes(regl, data.transform, data.positionCount)
            }),
            count: data.positionCount,
            instances: data.instanceCount,
            primitive: 'points'
        })
        return {
            draw: () => command(),
            get stats() {
                return command.stats
            },
            name: 'point'
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