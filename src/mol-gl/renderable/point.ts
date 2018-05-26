/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, DefineSpec, Values } from '../renderable/schema';
import { PointShaderCode } from '../shader-code';

export const PointSchema = {
    ...BaseSchema,
    aSize: AttributeSpec('float32', 1, 0),
    uSize: UniformSpec('f'),
    dSizeType: DefineSpec('string', ['uniform', 'attribute']),
    dPointSizeAttenuation: DefineSpec('boolean'),
}
export type PointSchema = typeof PointSchema
export type PointValues = Values<PointSchema>

export function PointRenderable(ctx: Context, values: PointValues, state: RenderableState): Renderable<PointValues> {
    const schema = { ...GlobalUniformSchema, ...PointSchema }
    const schaderCode = PointShaderCode
    const renderItem = createRenderItem(ctx, 'points', schaderCode, schema, values)

    return {
        draw: () => {
            renderItem.draw()
        },
        get values () { return values },
        get state () { return state },
        name: 'point',
        get program () { return renderItem.program },
        update: () => {
            renderItem.update()
        },
        dispose: () => {
            renderItem.destroy()
        }
    }
}