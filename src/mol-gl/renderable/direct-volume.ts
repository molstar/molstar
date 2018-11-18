/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { WebGLContext } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, GlobalUniformSchema, InternalSchema, TextureSpec, ValueSpec, ElementsSpec, DefineSpec, InternalValues } from './schema';
import { DirectVolumeShaderCode } from '../shader-code';
import { ValueCell } from 'mol-util';

export const DirectVolumeSchema = {
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
    boundingSphere: ValueSpec('sphere'),

    aPosition: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),

    uAlpha: UniformSpec('f'),
    uHighlightColor: UniformSpec('v3'),
    uSelectColor: UniformSpec('v3'),
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

    dGridTexType: DefineSpec('string', ['2d', '3d']),
    uGridTexDim: UniformSpec('v3'),
    tGridTex: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
}
export type DirectVolumeSchema = typeof DirectVolumeSchema
export type DirectVolumeValues = Values<DirectVolumeSchema>

export function DirectVolumeRenderable(ctx: WebGLContext, id: number, values: DirectVolumeValues, state: RenderableState): Renderable<DirectVolumeValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...DirectVolumeSchema }
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
        uPickable: ValueCell.create(state.pickable ? 1 : 0)
    }
    const shaderCode = DirectVolumeShaderCode
    const renderItem = createRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues })
    const renderable = createRenderable(renderItem, values, state);

    Object.defineProperty(renderable, 'opaque', { get: () => false });

    return renderable
}