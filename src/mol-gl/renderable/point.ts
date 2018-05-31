/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, DefineSpec, Values, InternalSchema } from '../renderable/schema';
import { PointShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

export const PointSchema = {
    ...BaseSchema,
    aSize: AttributeSpec('float32', 1, 0),
    uSize: UniformSpec('f'),
    dSizeType: DefineSpec('string', ['uniform', 'attribute']),
    dPointSizeAttenuation: DefineSpec('boolean'),
}
export type PointSchema = typeof PointSchema
export type PointValues = Values<PointSchema>

export function PointRenderable(ctx: Context, id: number, values: PointValues, state: RenderableState): Renderable<PointValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...PointSchema }
    const internalValues = {
        uObjectId: ValueCell.create(id)
    }
    const schaderCode = PointShaderCode
    const renderItem = createRenderItem(ctx, 'points', schaderCode, schema, { ...values, ...internalValues })

    return createRenderable(renderItem, values, state)
}