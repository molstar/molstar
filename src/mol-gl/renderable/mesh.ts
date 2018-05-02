/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { ColorData } from 'mol-geo/util/color-data';

import { Renderable } from '../renderable'
import { getBaseDefs, getBaseValues, getBaseDefines } from './util'
import { MeshShaderCode, addShaderDefines } from '../shader-code'
import { Context } from '../webgl/context';
import { createRenderItem, RenderItemProps, RenderItemState } from '../webgl/render-item';

type Mesh = 'mesh'

namespace Mesh {
    export type Props = {
        objectId: number
        alpha: number

        position: ValueCell<Float32Array>
        normal: ValueCell<Float32Array | undefined>
        id: ValueCell<Float32Array>

        color: ColorData
        transform: ValueCell<Float32Array>
        index: ValueCell<Uint32Array>

        indexCount: number
        instanceCount: number
        elementCount: number
        positionCount: number

        flatShaded?: boolean
        doubleSided?: boolean
    }

    export function create(ctx: Context, props: Props): Renderable<Props> {
        const defines = getBaseDefines(props)
        if (props.flatShaded) defines.FLAT_SHADED = ''
        if (props.doubleSided) defines.DOUBLE_SIDED = ''

        const defs: RenderItemProps = {
            ...getBaseDefs(props),
            shaderCode: addShaderDefines(defines, MeshShaderCode),
            drawMode: 'triangles',
            elementsKind: 'uint32'
        }
        const values: RenderItemState = {
            ...getBaseValues(props),
            drawCount: props.indexCount * 3,
            instanceCount: props.instanceCount,
            elements: props.index.ref.value
        }

        let renderItem = createRenderItem(ctx, defs, values)
        // let curProps = props

        return {
            draw: () => {
                renderItem.draw()
            },
            name: 'mesh',
            get program () { return renderItem.program },
            update: (newProps: Props) => {
                console.log('Updating mesh renderable')
            },
            dispose: () => {
                renderItem.destroy()
            }
        }
    }
}

export default Mesh