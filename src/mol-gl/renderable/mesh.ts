/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'

import { Renderable } from '../renderable'
import { ColorTexture } from '../util'
import { getBuffers, createTransformAttributes, fillSerial, createColorUniforms } from './util'
import Attribute from '../attribute';
import { MeshShaders } from '../shaders'

type Mesh = 'mesh'

type Uniforms = { [k: string]: REGL.Uniform | REGL.Texture }

namespace Mesh {
    export type Data = {
        position: ValueCell<Float32Array>
        normal: ValueCell<Float32Array>
        transform: ValueCell<Float32Array>
        color: ValueCell<ColorTexture>
        elements: ValueCell<Uint32Array>

        instanceCount: number
        elementCount: number
        positionCount: number
    }

    export function create(regl: REGL.Regl, data: Data, uniforms: Uniforms): Renderable {
        console.log(data)
        const instanceId = ValueCell.create(fillSerial(new Float32Array(data.instanceCount)))
        const command = regl({
            ...MeshShaders,
            uniforms: {
                objectId: uniforms.objectId || 0,
                instanceCount: data.instanceCount,
                ...createColorUniforms(regl, data.color),
                ...uniforms
            },
            attributes: getBuffers({
                instanceId: Attribute.create(regl, instanceId, data.instanceCount, { size: 1, divisor: 1 }),
                position: Attribute.create(regl, data.position, data.positionCount, { size: 3 }),
                normal: Attribute.create(regl, data.normal, data.positionCount, { size: 3 }),
                ...createTransformAttributes(regl, data.transform, data.instanceCount)
            }),
            elements: regl.elements({
                data: data.elements.ref.value,
                primitive: 'triangles',
                type: 'uint32',
                count: data.elementCount * 3,
                // length: count * 3 * 2
            }),
            instances: data.instanceCount,
        })
        return {
            draw: () => {
                command()
                console.log(command.stats)
            }
        }
    }
}

export default Mesh