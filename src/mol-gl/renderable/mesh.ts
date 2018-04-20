/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'
import { ColorData } from 'mol-geo/color';

import { Renderable } from '../renderable'
import { createBaseDefines, createBaseUniforms, createBaseAttributes, destroyAttributes, destroyUniforms } from './util'
import { MeshShaders, addDefines } from '../shaders'

type Mesh = 'mesh'

namespace Mesh {
    export type Data = {
        objectId: number

        position: ValueCell<Float32Array>
        normal?: ValueCell<Float32Array>
        id: ValueCell<Float32Array>

        color: ColorData
        transform: ValueCell<Float32Array>
        index: ValueCell<Uint32Array>

        indexCount: number
        instanceCount: number
        elementCount: number
        positionCount: number
    }

    export function create(regl: REGL.Regl, props: Data): Renderable {
        const defines = createBaseDefines(regl, props)
        const uniforms = createBaseUniforms(regl, props)
        const attributes = createBaseAttributes(regl, props)
        const elements = regl.elements({
            data: props.index.ref.value,
            primitive: 'triangles',
            type: 'uint32',
            count: props.indexCount * 3
        })

        const command = regl({
            ...addDefines(defines, MeshShaders),
            uniforms,
            attributes,
            elements,
            instances: props.instanceCount,
        })
        return {
            draw: () => {
                command()
            },
            get stats() {
                return command.stats
            },
            name: 'mesh',
            dispose: () => {
                destroyAttributes(attributes)
                destroyUniforms(uniforms)
                elements.destroy()
            }
        }
    }
}

export default Mesh