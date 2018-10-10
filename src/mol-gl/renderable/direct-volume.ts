/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, GlobalUniformSchema, InternalSchema, TextureSpec, ValueSpec, ElementsSpec, DefineSpec, InternalValues } from './schema';
import { DirectVolumeShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

const DirectVolumeBaseSchema = {
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
    dRenderMode: DefineSpec('string', ['isosurface', 'volume']),
    tTransferTex: TextureSpec('image-uint8', 'rgba', 'ubyte', 'linear'),
}

function getInternalValues(ctx: Context, id: number): InternalValues {
    return {
        dWebGL2: ValueCell.create(ctx.isWebGL2),
        uObjectId: ValueCell.create(id)
    }
}

//

export const DirectVolume2dSchema = {
    ...DirectVolumeBaseSchema,
    uGridTexDim: UniformSpec('v2'),
    tGridTex: TextureSpec('image-uint8', 'rgba', 'ubyte', 'linear'),
}
export type DirectVolume2dSchema = typeof DirectVolume2dSchema
export type DirectVolume2dValues = Values<DirectVolume2dSchema>

export function DirectVolume2dRenderable(ctx: Context, id: number, values: DirectVolume2dValues, state: RenderableState): Renderable<DirectVolume2dValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...DirectVolume2dSchema }
    const internalValues = getInternalValues(ctx, id)
    const shaderCode = DirectVolumeShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues })
    const renderable = createRenderable(renderItem, values, state);

    Object.defineProperty(renderable, 'opaque', { get: () => false });

    return renderable
}