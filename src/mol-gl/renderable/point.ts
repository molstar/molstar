/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { Renderable, AttributesMutator, AttributesData } from '../renderable'
import { createAttributes, getBuffers, getData } from './util'

const pointVert = require('mol-gl/shader/point.vert')
const pointFrag = require('mol-gl/shader/point.frag')

type Point = 'point'

namespace Point {
    export interface Data extends AttributesData {
        position: Float32Array
    }
    export function create(regl: REGL.Regl, data: Data): Renderable<Data> {
        const attributes = createAttributes(regl, data)
        const command = regl({
            vert: pointVert,
            frag: pointFrag,
            attributes: getBuffers(attributes),
            count: data.position.length / 3,
            primitive: 'points'
        })
        return {
            draw: () => command(),
            update: (mutator: AttributesMutator<Data>) => {
                mutator(getData(attributes))
                for (const k of Object.keys(attributes)) {
                    attributes[k].reload()
                }
            }
        }
    }
}

export default Point