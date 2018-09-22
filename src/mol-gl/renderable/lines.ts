/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, DefineSpec, Values, InternalSchema, SizeSchema, ElementsSpec } from './schema';
import { ValueCell } from 'mol-util';
import { LinesShaderCode } from '../shader-code';

export const LinesSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aMapping: AttributeSpec('float32', 2, 0),
    aStart: AttributeSpec('float32', 3, 0),
    aEnd: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),
    dLineSizeAttenuation: DefineSpec('boolean'),
    dDoubleSided: DefineSpec('boolean'),
    dFlipSided: DefineSpec('boolean'),
}
export type LinesSchema = typeof LinesSchema
export type LinesValues = Values<LinesSchema>

export function LinesRenderable(ctx: Context, id: number, values: LinesValues, state: RenderableState): Renderable<LinesValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...LinesSchema }
    const internalValues = {
        uObjectId: ValueCell.create(id)
    }
    const shaderCode = LinesShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues })

    return createRenderable(renderItem, values, state);
}