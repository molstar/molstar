/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { Renderable } from '../renderable'
import { getBaseValues, getBaseDefs, getBaseDefines } from './util'
import { PointShaderCode, addShaderDefines } from '../shader-code'
import { ColorData } from 'mol-geo/util/color-data';
import { SizeData } from 'mol-geo/util/size-data';
import { Context } from '../webgl/context';
import { createRenderItem, RenderItemState, RenderItemProps } from '../webgl/render-item';

type Point = 'point'

namespace Point {
    export type Props = {
        objectId: number

        position: ValueCell<Float32Array>
        id: ValueCell<Float32Array>

        size: SizeData
        color: ColorData
        transform: ValueCell<Float32Array>

        instanceCount: number
        elementCount: number
        positionCount: number,

        usePointSizeAttenuation?: boolean
    }

    export function create<T = Props>(ctx: Context, props: Props): Renderable<Props> {
        const defines = getBaseDefines(props)
        if (props.usePointSizeAttenuation) defines.POINT_SIZE_ATTENUATION = ''

        const defs: RenderItemProps = {
            ...getBaseDefs(props),
            shaderCode: addShaderDefines(defines, PointShaderCode),
            drawMode: 'points'
        }
        const values: RenderItemState = {
            ...getBaseValues(props),
            drawCount: props.positionCount,
            instanceCount: props.instanceCount
        }

        let renderItem = createRenderItem(ctx, defs, values)
        // let curProps = props

        return {
            draw: () => {
                renderItem.draw()
            },
            name: 'point',
            get program () { return renderItem.program },
            update: (newProps: Props) => {
                console.log('Updating point renderable')
            },
            dispose: () => {
                renderItem.destroy()
            }
        }
    }
}

export default Point