/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, DefineSpec, Values, InternalSchema, InternalValues, UniformSpec, TextureSpec } from './schema';
import { MeshShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

export const IsosurfaceSchema = {
    ...BaseSchema,

    aIndex: AttributeSpec('float32', 1, 0),
    aNormal: AttributeSpec('float32', 3, 0),

    uPositionTexDim: UniformSpec('v2'),
    tPosition: TextureSpec('texture', 'rgba', 'float', 'nearest'),

    dFlatShaded: DefineSpec('boolean'),
    dDoubleSided: DefineSpec('boolean'),
    dFlipSided: DefineSpec('boolean'),
    dPositionTexture: DefineSpec('boolean'),
}
export type IsosurfaceSchema = typeof IsosurfaceSchema
export type IsosurfaceValues = Values<IsosurfaceSchema>

export function IsosurfaceRenderable(ctx: WebGLContext, id: number, values: IsosurfaceValues, state: RenderableState, materialId: number): Renderable<IsosurfaceValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...IsosurfaceSchema }
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
        uPickable: ValueCell.create(state.pickable ? 1 : 0)
    }
    const shaderCode = MeshShaderCode
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId)

    return createRenderable(renderItem, values, state)
}