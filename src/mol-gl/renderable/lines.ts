/**
 * Copyright (c) 2018-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem, Transparency } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, DefineSpec, Values, InternalSchema, SizeSchema, ElementsSpec, InternalValues, GlobalTextureSchema, UniformSpec } from './schema';
import { ValueCell } from '../../mol-util';
import { LinesShaderCode } from '../shader-code';

export const LinesSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aGroup: AttributeSpec('float32', 1, 0),
    aMapping: AttributeSpec('float32', 2, 0),
    aStart: AttributeSpec('float32', 3, 0),
    aEnd: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),
    dLineSizeAttenuation: DefineSpec('boolean'),
    uDoubleSided: UniformSpec('b', 'material'),
    dFlipSided: DefineSpec('boolean'),
};
export type LinesSchema = typeof LinesSchema
export type LinesValues = Values<LinesSchema>

export function LinesRenderable(ctx: WebGLContext, id: number, values: LinesValues, state: RenderableState, materialId: number, transparency: Transparency): Renderable<LinesValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...LinesSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = LinesShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId, transparency);

    return createRenderable(renderItem, values, state);
}