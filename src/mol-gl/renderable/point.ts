/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, DefineSpec, Values, InternalSchema, TextureSpec } from './schema';
import { PointShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

export const PointSchema = {
    ...BaseSchema,
    aSize: AttributeSpec('float32', 1, 0),
    uSize: UniformSpec('f'),
    uSizeTexDim: UniformSpec('v2'),
    tSize: TextureSpec('alpha', 'ubyte'),
    dSizeType: DefineSpec('string', ['uniform', 'attribute']),
    dPointSizeAttenuation: DefineSpec('boolean'),
    dPointFilledCircle: DefineSpec('boolean'),
    uPointEdgeBleach: UniformSpec('f'),
}
export type PointSchema = typeof PointSchema
export type PointValues = Values<PointSchema>

export function PointRenderable(ctx: Context, id: number, values: PointValues, state: RenderableState): Renderable<PointValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...PointSchema }
    const internalValues = {
        uObjectId: ValueCell.create(id)
    }
    const shaderCode = PointShaderCode
    const renderItem = createRenderItem(ctx, 'points', shaderCode, schema, { ...values, ...internalValues })
    const renderable = createRenderable(renderItem, values, state);

    const isOpaque = Object.getOwnPropertyDescriptor(renderable, 'opaque')!.get as () => boolean
    Object.defineProperty(renderable, 'opaque', {
        get: () => isOpaque() && !values.dPointFilledCircle.ref.value && values.uPointEdgeBleach.ref.value === 0
    });

    return renderable
}