/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem, GraphicsRenderVariant } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, Values, InternalSchema, SizeSchema, InternalValues, ElementsSpec, ValueSpec, DefineSpec, GlobalTextureSchema, UniformSpec } from './schema';
import { CylindersShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const CylindersSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aGroup: AttributeSpec('float32', 1, 0),
    aStart: AttributeSpec('float32', 3, 0),
    aEnd: AttributeSpec('float32', 3, 0),
    aMapping: AttributeSpec('float32', 3, 0),
    aScale: AttributeSpec('float32', 1, 0),
    aCap: AttributeSpec('float32', 1, 0),
    elements: ElementsSpec('uint32'),

    padding: ValueSpec('number'),
    uDoubleSided: UniformSpec('b', 'material'),
    dIgnoreLight: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('boolean'),
    dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
    dSolidInterior: DefineSpec('boolean'),
    uBumpFrequency: UniformSpec('f', 'material'),
    uBumpAmplitude: UniformSpec('f', 'material'),
};
export type CylindersSchema = typeof CylindersSchema
export type CylindersValues = Values<CylindersSchema>

export function CylindersRenderable(ctx: WebGLContext, id: number, values: CylindersValues, state: RenderableState, materialId: number, variants: GraphicsRenderVariant[]): Renderable<CylindersValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...CylindersSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = CylindersShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId, variants);
    return createRenderable(renderItem, values, state);
}