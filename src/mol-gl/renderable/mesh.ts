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

// TODO
interface Elements {

}

namespace Mesh {
    export type DataType = {
        position: { type: Float32Array, itemSize: 3 }
        offset: { type: Float32Array, itemSize: 3 }
    }
    export type Data = { [K in keyof DataType]: DataType[K]['type'] }
    export type Attributes = { [K in keyof Data]: Attribute<Data[K]> }

    export function create(regl: REGL.Regl, attributes: Attributes, elements?: Elements): Renderable<Data> {
        console.log('mesh', {
            count: attributes.position.getCount(),
            instances: attributes.transform.getCount(),
        })
        const command = regl({
            ...MeshShaders,
            attributes: getBuffers(attributes),
            count: attributes.position.getCount(),
            instances: attributes.transform.getCount(),
            primitive: 'triangles'
        })
        return {
            draw: () => command(),
        }
    }
}

export default Mesh