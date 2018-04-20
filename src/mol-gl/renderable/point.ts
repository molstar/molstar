/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import REGL = require('regl');
import { ValueCell } from 'mol-util/value-cell'

import { Renderable } from '../renderable'
import { createBaseDefines, createBaseUniforms, createBaseAttributes, destroyUniforms, destroyAttributes } from './util'
import { PointShaders, addDefines } from '../shaders'
import { ColorData } from 'mol-geo/color';

type Point = 'point'

namespace Point {
    export type Data = {
        objectId: number

        position: ValueCell<Float32Array>
        size?: ValueCell<Float32Array>
        id: ValueCell<Float32Array>

        color: ColorData
        transform: ValueCell<Float32Array>

        instanceCount: number
        elementCount: number
        positionCount: number
    }

    export function create(regl: REGL.Regl, props: Data): Renderable {
        const defines = createBaseDefines(regl, props)
        const uniforms = createBaseUniforms(regl, props)
        const attributes = createBaseAttributes(regl, props)

        const command = regl({
            ...addDefines(defines, PointShaders),
            uniforms,
            attributes,
            count: props.positionCount,
            instances: props.instanceCount,
            primitive: 'points'
        })
        return {
            draw: () => command(),
            get stats() {
                return command.stats
            },
            name: 'point',
            dispose: () => {
                destroyAttributes(attributes)
                destroyUniforms(uniforms)
            }
        }
    }
}

export default Point