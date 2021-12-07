/**
 * Copyright (c) 2019-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, DefineSpec, Values, InternalSchema, InternalValues, UniformSpec, TextureSpec, GlobalTextureSchema, ValueSpec } from './schema';
import { MeshShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const TextureMeshSchema = {
    ...BaseSchema,
    uGeoTexDim: UniformSpec('v2', 'buffered'),
    tPosition: TextureSpec('texture', 'rgb', 'float', 'nearest'),
    tGroup: TextureSpec('texture', 'alpha', 'float', 'nearest'),
    tNormal: TextureSpec('texture', 'rgb', 'float', 'nearest'),

    dFlatShaded: DefineSpec('boolean'),
    dDoubleSided: DefineSpec('boolean'),
    dFlipSided: DefineSpec('boolean'),
    dIgnoreLight: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('boolean'),
    dGeoTexture: DefineSpec('boolean'),
    uBumpFrequency: UniformSpec('f'),
    uBumpAmplitude: UniformSpec('f'),
    meta: ValueSpec('unknown')
};
export type TextureMeshSchema = typeof TextureMeshSchema
export type TextureMeshValues = Values<TextureMeshSchema>

export function TextureMeshRenderable(ctx: WebGLContext, id: number, values: TextureMeshValues, state: RenderableState, materialId: number): Renderable<TextureMeshValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...TextureMeshSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = MeshShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId);

    return createRenderable(renderItem, values, state);
}