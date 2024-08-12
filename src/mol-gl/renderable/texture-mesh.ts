/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem, Transparency } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, DefineSpec, Values, InternalSchema, InternalValues, UniformSpec, TextureSpec, GlobalTextureSchema, ValueSpec } from './schema';
import { MeshShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const TextureMeshSchema = {
    ...BaseSchema,
    uGeoTexDim: UniformSpec('v2', 'buffered'),
    tPosition: TextureSpec('texture', 'rgb', 'float', 'nearest'),
    tGroup: TextureSpec('texture', 'alpha', 'float', 'nearest'),
    tNormal: TextureSpec('texture', 'rgb', 'float', 'nearest'),
    dVaryingGroup: DefineSpec('boolean'),
    dFlatShaded: DefineSpec('boolean'),
    uDoubleSided: UniformSpec('b', 'material'),
    dFlipSided: DefineSpec('boolean'),
    dIgnoreLight: DefineSpec('boolean'),
    dCelShaded: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('string', ['off', 'on', 'inverted']),
    dTransparentBackfaces: DefineSpec('string', ['off', 'on', 'opaque']),
    uBumpFrequency: UniformSpec('f', 'material'),
    uBumpAmplitude: UniformSpec('f', 'material'),
    meta: ValueSpec('unknown')
};
export type TextureMeshSchema = typeof TextureMeshSchema
export type TextureMeshValues = Values<TextureMeshSchema>

export function TextureMeshRenderable(ctx: WebGLContext, id: number, values: TextureMeshValues, state: RenderableState, materialId: number, transparency: Transparency): Renderable<TextureMeshValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...TextureMeshSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = MeshShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId, transparency);

    return createRenderable(renderItem, values, state);
}