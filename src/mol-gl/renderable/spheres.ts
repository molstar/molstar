/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem, Transparency } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, Values, InternalSchema, SizeSchema, InternalValues, ValueSpec, DefineSpec, GlobalTextureSchema, UniformSpec, TextureSpec } from './schema';
import { SpheresShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const SpheresSchema = {
    ...BaseSchema,
    ...SizeSchema,

    uTexDim: UniformSpec('v2'),
    tPositionGroup: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),

    padding: ValueSpec('number'),
    uDoubleSided: UniformSpec('b', 'material'),
    dIgnoreLight: DefineSpec('boolean'),
    dCelShaded: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
    dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
    dSolidInterior: DefineSpec('boolean'),
    dClipPrimitive: DefineSpec('boolean'),
    dApproximate: DefineSpec('boolean'),
    uAlphaThickness: UniformSpec('f'),
    uBumpFrequency: UniformSpec('f', 'material'),
    uBumpAmplitude: UniformSpec('f', 'material'),

    lodLevels: ValueSpec('unknown'),
    centerBuffer: ValueSpec('float32'),
    groupBuffer: ValueSpec('float32'),
};
export type SpheresSchema = typeof SpheresSchema
export type SpheresValues = Values<SpheresSchema>

export function SpheresRenderable(ctx: WebGLContext, id: number, values: SpheresValues, state: RenderableState, materialId: number, transparency: Transparency): Renderable<SpheresValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...SpheresSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = SpheresShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId, transparency);
    return createRenderable(renderItem, values, state);
}