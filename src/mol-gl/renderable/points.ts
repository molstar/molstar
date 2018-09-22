/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, DefineSpec, Values, InternalSchema, SizeSchema } from './schema';
import { PointsShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

export const PointsSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aPosition: AttributeSpec('float32', 3, 0),
    dPointSizeAttenuation: DefineSpec('boolean'),
    dPointFilledCircle: DefineSpec('boolean'),
    uPointEdgeBleach: UniformSpec('f'),
}
export type PointsSchema = typeof PointsSchema
export type PointsValues = Values<PointsSchema>

export function PointsRenderable(ctx: Context, id: number, values: PointsValues, state: RenderableState): Renderable<PointsValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...PointsSchema }
    const internalValues = {
        uObjectId: ValueCell.create(id)
    }
    const shaderCode = PointsShaderCode
    const renderItem = createRenderItem(ctx, 'points', shaderCode, schema, { ...values, ...internalValues })
    const renderable = createRenderable(renderItem, values, state);

    const isOpaque = Object.getOwnPropertyDescriptor(renderable, 'opaque')!.get as () => boolean
    Object.defineProperty(renderable, 'opaque', {
        get: () => isOpaque() && !values.dPointFilledCircle.ref.value && values.uPointEdgeBleach.ref.value === 0
    });

    return renderable
}