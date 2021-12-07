/**
 * Copyright (c) 2020-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
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
    dDoubleSided: DefineSpec('boolean'),
    dIgnoreLight: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('boolean'),
    uBumpFrequency: UniformSpec('f'),
    uBumpAmplitude: UniformSpec('f'),
};
export type CylindersSchema = typeof CylindersSchema
export type CylindersValues = Values<CylindersSchema>

export function CylindersRenderable(ctx: WebGLContext, id: number, values: CylindersValues, state: RenderableState, materialId: number): Renderable<CylindersValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...CylindersSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = CylindersShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId);
    return createRenderable(renderItem, values, state);
}