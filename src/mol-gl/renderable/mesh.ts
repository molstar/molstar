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

type AttributeColor = { type: 'attribute', value: ValueCell<Float32Array> }
type InstanceColor = { type: 'instance', value: ValueCell<ColorTexture> }
type ElementColor = { type: 'element', value: ValueCell<ColorTexture> }
export type Color = AttributeColor | InstanceColor | ElementColor

namespace Mesh {
    export type Data = {
        position: ValueCell<Float32Array>
        normal: ValueCell<Float32Array>
        readonly color: Color
        transform: ValueCell<Float32Array>
        elements: ValueCell<Uint32Array>

        instanceCount: number
        elementCount: number
        positionCount: number
    }

    export function create(regl: REGL.Regl, data: Data, _uniforms: Uniforms): Renderable {
        const instanceId = ValueCell.create(fillSerial(new Float32Array(data.instanceCount)))
        const uniforms = {
            objectId: _uniforms.objectId || 0,
            instanceCount: data.instanceCount,
            ..._uniforms
        }
        if (data.color.type === 'instance' || data.color.type === 'element') {
            Object.assign(uniforms, createColorUniforms(regl, data.color.value as ValueCell<ColorTexture>))
        }
        const attributes = getBuffers({
            instanceId: Attribute.create(regl, instanceId, data.instanceCount, { size: 1, divisor: 1 }),
            position: Attribute.create(regl, data.position, data.positionCount, { size: 3 }),
            normal: Attribute.create(regl, data.normal, data.positionCount, { size: 3 }),
            ...createTransformAttributes(regl, data.transform, data.instanceCount)
        })
        if (data.color.type === 'attribute') {
            attributes.color = Attribute.create(regl, data.color.value as ValueCell<Float32Array>, data.positionCount, { size: 3 }).buffer
        }
        const command = regl({
            ...MeshShaders,
            uniforms,
            attributes,
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
            },
            get stats() {
                return command.stats
            },
            name: 'mesh'
        }
    }
}

export default Mesh