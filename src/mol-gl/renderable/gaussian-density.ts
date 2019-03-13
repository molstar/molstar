/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { WebGLContext } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, ValueSpec, DefineSpec, TextureSpec } from './schema';
import { GaussianDensityShaderCode } from '../shader-code';

export const GaussianDensitySchema = {
    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    aRadius: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),
    aGroup: AttributeSpec('float32', 1, 0),

    uCurrentSlice: UniformSpec('f'),
    uCurrentX: UniformSpec('f'),
    uCurrentY: UniformSpec('f'),
    uBboxMin: UniformSpec('v3'),
    uBboxMax: UniformSpec('v3'),
    uBboxSize: UniformSpec('v3'),
    uGridDim: UniformSpec('v3'),
    uGridTexDim: UniformSpec('v3'),
    uAlpha: UniformSpec('f'),
    tMinDistanceTex: TextureSpec('texture', 'rgba', 'ubyte', 'nearest'),

    dGridTexType: DefineSpec('string', ['2d', '3d']),
    dCalcType: DefineSpec('string', ['density', 'minDistance', 'groupId']),
}
export type GaussianDensitySchema = typeof GaussianDensitySchema
export type GaussianDensityValues = Values<GaussianDensitySchema>

export function GaussianDensityRenderable(ctx: WebGLContext, id: number, values: GaussianDensityValues, state: RenderableState): Renderable<GaussianDensityValues> {
    const schema = { ...GaussianDensitySchema }
    const shaderCode = GaussianDensityShaderCode
    const renderItem = createRenderItem(ctx, 'points', shaderCode, schema, values, -1)

    return createRenderable(renderItem, values, state);
}