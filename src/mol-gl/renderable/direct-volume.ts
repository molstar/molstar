/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, GlobalUniformSchema, InternalSchema, TextureSpec, ValueSpec, ElementsSpec, DefineSpec } from './schema';
import { DirectVolumeShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

export const DirectVolumeSchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    aPosition: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),

    uAlpha: UniformSpec('f'),
    dUseFog: DefineSpec('boolean'),

    uIsoValue: UniformSpec('f'),
    uBboxMin: UniformSpec('v3'),
    uBboxMax: UniformSpec('v3'),
    uBboxSize: UniformSpec('v3'),
    uTransform: UniformSpec('m4'),
    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v2'),
    tGridTex: TextureSpec('rgba', 'ubyte', 'linear'),
    dRenderMode: DefineSpec('string', ['isosurface', 'volume']),
    tTransferTex: TextureSpec('rgba', 'ubyte', 'linear'),
}
export type DirectVolumeSchema = typeof DirectVolumeSchema
export type DirectVolumeValues = Values<DirectVolumeSchema>

export function DirectVolumeRenderable(ctx: Context, id: number, values: DirectVolumeValues, state: RenderableState): Renderable<DirectVolumeValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...DirectVolumeSchema }
    const internalValues = {
        uObjectId: ValueCell.create(id)
    }
    const shaderCode = DirectVolumeShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues })
    const renderable = createRenderable(renderItem, values, state);

    Object.defineProperty(renderable, 'opaque', { get: () => false });

    return renderable
}