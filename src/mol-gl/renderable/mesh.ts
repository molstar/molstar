/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'
import { ColorData } from 'mol-geo/color';

import { Renderable } from '../renderable'
import { getBuffers, createTransformAttributes, fillSerial, createColorUniforms } from './util'
import Attribute from '../attribute';
import { MeshShaders, addDefines, ShaderDefines } from '../shaders'

type Mesh = 'mesh'

type Uniforms = { [k: string]: REGL.Uniform | REGL.Texture }

function getColorDefines(color: ColorData) {
    const defines: ShaderDefines = {}
    switch (color.type) {
        case 'uniform': defines.UNIFORM_COLOR = ''; break;
        case 'attribute': defines.ATTRIBUTE_COLOR = ''; break;
        case 'element': defines.ELEMENT_COLOR = ''; break;
        case 'instance': defines.INSTANCE_COLOR = ''; break;
        case 'element-instance': defines.ELEMENT_INSTANCE_COLOR = ''; break;
    }
    return defines
}

namespace Mesh {
    export type Data = {
        position: ValueCell<Float32Array>
        normal: ValueCell<Float32Array>
        id: ValueCell<Float32Array>

        readonly color: ColorData
        transform: ValueCell<Float32Array>
        index: ValueCell<Uint32Array>

        indexCount: number
        instanceCount: number
        elementCount: number
        positionCount: number
    }

    export function create(regl: REGL.Regl, data: Data, _uniforms: Uniforms): Renderable {
        const defines = getColorDefines(data.color)
        const instanceId = ValueCell.create(fillSerial(new Float32Array(data.instanceCount)))
        const uniforms = {
            objectId: _uniforms.objectId || 0,
            instanceCount: data.instanceCount,
            elementCount: data.elementCount,
            ..._uniforms
        }
        if (data.color.type === 'instance' || data.color.type === 'element' || data.color.type === 'element-instance') {
            Object.assign(uniforms, createColorUniforms(regl, data.color.value))
        } else if (data.color.type === 'uniform') {
            Object.assign(uniforms, { color: data.color.value })
        }
        const attributes = getBuffers({
            instanceId: Attribute.create(regl, instanceId, data.instanceCount, { size: 1, divisor: 1 }),
            position: Attribute.create(regl, data.position, data.positionCount, { size: 3 }),
            normal: Attribute.create(regl, data.normal, data.positionCount, { size: 3 }),

            elementId: Attribute.create(regl, data.id, data.positionCount, { size: 1 }),
            ...createTransformAttributes(regl, data.transform, data.instanceCount)
        })
        if (data.color.type === 'attribute') {
            attributes.color = Attribute.create(regl, data.color.value, data.positionCount, { size: 3 }).buffer
        }
        const command = regl({
            ...addDefines(defines, MeshShaders),
            uniforms,
            attributes,
            elements: regl.elements({
                data: data.index.ref.value,
                primitive: 'triangles',
                type: 'uint32',
                count: data.indexCount * 3
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