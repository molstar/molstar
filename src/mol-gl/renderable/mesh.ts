/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { ColorData } from 'mol-geo/util/color-data';

import { Renderable, BaseProps } from '../renderable'
import { getBaseDefs, getBaseValues, getBaseDefines, updateBaseValues } from './util'
import { MeshShaderCode, addShaderDefines } from '../shader-code'
import { Context } from '../webgl/context';
import { createRenderItem, RenderItemProps, RenderItemState } from '../webgl/render-item';
import { deepEqual } from 'mol-util';

type Mesh = 'mesh'

namespace Mesh {
    export type Props = {
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
    } & BaseProps

    function getDefs(props: Props) {
        const defines = getBaseDefines(props)
        if (props.flatShaded) defines.FLAT_SHADED = ''
        if (props.doubleSided) defines.DOUBLE_SIDED = ''
        if (props.flipSided) defines.FLIP_SIDED = ''

        const defs: RenderItemProps = {
            ...getBaseDefs(props),
            shaderCode: addShaderDefines(defines, MeshShaderCode),
            drawMode: 'triangles',
            elementsKind: 'uint32'
        }
        return defs
    }

    function getVals(props: Props) {
        const vals: RenderItemState = {
            ...getBaseValues(props),
            drawCount: ValueCell.create(props.indexCount * 3),
            instanceCount: ValueCell.create(props.instanceCount),
            elements: props.index.ref.value
        }
        return vals
    }

    function updateVals(vals: RenderItemState, props: Props) {
        updateBaseValues(vals, props)
        if (props.instanceCount !== vals.instanceCount.ref.value) {
            ValueCell.update(vals.instanceCount, props.instanceCount)
        }
        const drawCount = props.indexCount * 3
        if (drawCount !== vals.drawCount.ref.value) {
            ValueCell.update(vals.drawCount, drawCount)
        }
    }

    export function create(ctx: Context, props: Props): Renderable<Props> {
        let curDefs = getDefs(props)
        let curVals = getVals(props)
        let renderItem = createRenderItem(ctx, curDefs, curVals)

        return {
            draw: () => {
                renderItem.draw()
            },
            name: 'mesh',
            get program () { return renderItem.program },
            update: (newProps: Props) => {
                const newDefs = getDefs(props)
                if (deepEqual(curDefs, newDefs)) {
                    updateVals(curVals, props)
                    renderItem.update()
                } else {
                    console.log('mesh defs changed, destroy and rebuild render-item')
                    renderItem.destroy()
                    curVals = getVals(props)
                    curDefs = newDefs
                    renderItem = createRenderItem(ctx, curDefs, curVals)
                }
            },
            dispose: () => {
                renderItem.destroy()
            }
        }
    }
}

export default Mesh