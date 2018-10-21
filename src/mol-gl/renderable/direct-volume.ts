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

export const DirectVolumeBaseSchema = {
    aColor: AttributeSpec('float32', 3, 0), // TODO not used, just for type checking
    uColor: UniformSpec('v3'),
    uColorTexDim: UniformSpec('v2'),
    tColor: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    dColorType: DefineSpec('string', ['uniform', 'instance', 'group', 'group_instance']),

    uMarkerTexDim: UniformSpec('v2'),
    tMarker: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),

    uInstanceCount: UniformSpec('i'),
    uGroupCount: UniformSpec('i'),

    aInstance: AttributeSpec('float32', 1, 1),
    aTransform: AttributeSpec('float32', 16, 1),

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
    dMaxSteps: DefineSpec('number'),
    uTransform: UniformSpec('m4'),
    uGridDim: UniformSpec('v3'),
    dRenderMode: DefineSpec('string', ['isosurface', 'volume']),
    tTransferTex: TextureSpec('image-uint8', 'rgba', 'ubyte', 'linear'),
}
export type DirectVolumeBaseSchema = typeof DirectVolumeBaseSchema
export type DirectVolumeBaseValues = Values<DirectVolumeBaseSchema>

function getInternalValues(ctx: Context, id: number): InternalValues {
    return {
        uObjectId: ValueCell.create(id)
    }
}

function DirectVolumeRenderable<T extends DirectVolumeBaseValues, S extends DirectVolumeBaseSchema>(ctx: Context, id: number, values: T, state: RenderableState, schema: S): Renderable<T> {
    const fullSchema = Object.assign({}, GlobalUniformSchema, InternalSchema, schema)
    const internalValues = getInternalValues(ctx, id)
    const fullValues = Object.assign({}, values, internalValues)
    const shaderCode = DirectVolumeShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', shaderCode, fullSchema, fullValues)
    const renderable = createRenderable(renderItem, values, state);

    Object.defineProperty(renderable, 'opaque', { get: () => true });

    return renderable
}

// via 2d texture

export const DirectVolume2dSchema = {
    ...DirectVolumeBaseSchema,
    dGridTexType: DefineSpec('string', ['2d']),
    uGridTexDim: UniformSpec('v2'),
    tGridTex: TextureSpec('texture2d', 'rgba', 'ubyte', 'linear'),
}
export type DirectVolume2dSchema = typeof DirectVolume2dSchema
export type DirectVolume2dValues = Values<DirectVolume2dSchema>

export function DirectVolume2dRenderable(ctx: Context, id: number, values: DirectVolume2dValues, state: RenderableState): Renderable<DirectVolume2dValues> {
    return DirectVolumeRenderable(ctx, id, values, state, DirectVolume2dSchema)
}

// via 3d texture

export const DirectVolume3dSchema = {
    ...DirectVolumeBaseSchema,
    dGridTexType: DefineSpec('string', ['3d']),
    tGridTex: TextureSpec('texture3d', 'rgba', 'ubyte', 'linear'),
}
export type DirectVolume3dSchema = typeof DirectVolume3dSchema
export type DirectVolume3dValues = Values<DirectVolume3dSchema>

export function DirectVolume3dRenderable(ctx: Context, id: number, values: DirectVolume3dValues, state: RenderableState): Renderable<DirectVolume3dValues> {
    return DirectVolumeRenderable(ctx, id, values, state, DirectVolume3dSchema)
}