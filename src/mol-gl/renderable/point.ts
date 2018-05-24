/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'

import { Renderable, BaseProps } from '../renderable'
import { getBaseValues, getBaseDefs, getBaseDefines } from './util'
import { PointShaderCode, addShaderDefines } from '../shader-code'
import { ColorData } from 'mol-geo/util/color-data';
import { SizeData } from 'mol-geo/util/size-data';
import { Context } from '../webgl/context';
import { createRenderItem, RenderItemState, RenderItemProps } from '../webgl/render-item';

type Point = 'point'

namespace Point {
    export type Props = {
        position: ValueCell<Float32Array>
        id: ValueCell<Float32Array>

        size: SizeData
        color: ColorData
        transform: ValueCell<Float32Array>

        instanceCount: number
        elementCount: number
        positionCount: number,

        usePointSizeAttenuation?: boolean
    } & BaseProps

    function getDefs(props: Props) {
        const defines = getBaseDefines(props)
        if (props.usePointSizeAttenuation) defines.POINT_SIZE_ATTENUATION = ''

        const defs: RenderItemProps = {
            ...getBaseDefs(props),
            shaderCode: addShaderDefines(defines, PointShaderCode),
            drawMode: 'points'
        }
        return defs
    }

    function getVals(props: Props) {
        const vals: RenderItemState = {
            ...getBaseValues(props),
            drawCount: ValueCell.create(props.positionCount),
            instanceCount: ValueCell.create(props.instanceCount)
        }
        return vals
    }

    function getRenderItem(ctx: Context, props: Props) {
        const defs = getDefs(props)
        const vals = getVals(props)
        return createRenderItem(ctx, defs, vals)
    }

    export function create<T = Props>(ctx: Context, props: Props): Renderable<Props> {
        // const defs = getDefs(props)

        let renderItem = getRenderItem(ctx, props)
        // let curProps = props

        return {
            draw: () => {
                renderItem.draw()
            },
            name: 'point',
            get program () { return renderItem.program },
            update: (newProps: Props) => {
                console.log('Updating point renderable')
                renderItem.destroy()
                renderItem = getRenderItem(ctx, { ...props, ...newProps })
            },
            dispose: () => {
                renderItem.destroy()
            }
        }
    }
}

export default Point