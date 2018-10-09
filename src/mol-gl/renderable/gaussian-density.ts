/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable'
import { Context } from '../webgl/context';
import { createRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, ValueSpec, DefineSpec } from './schema';
import { GaussianDensityShaderCode } from '../shader-code';

export const GaussianDensitySchema = {
    dWebGL2: DefineSpec('boolean'),

    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    aRadius: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),

    uCurrentSlice: UniformSpec('f'),
    uCurrentX: UniformSpec('f'),
    uCurrentY: UniformSpec('f'),
    uBboxMin: UniformSpec('v3'),
    uBboxMax: UniformSpec('v3'),
    uBboxSize: UniformSpec('v3'),
    uGridDim: UniformSpec('v3'),
    uAlpha: UniformSpec('f'),
}
export type GaussianDensitySchema = typeof GaussianDensitySchema
export type GaussianDensityValues = Values<GaussianDensitySchema>

export function GaussianDensityRenderable(ctx: Context, id: number, values: GaussianDensityValues, state: RenderableState): Renderable<GaussianDensityValues> {
    const schema = { ...GaussianDensitySchema }
    const shaderCode = GaussianDensityShaderCode
    const renderItem = createRenderItem(ctx, 'points', shaderCode, schema, values)

    return createRenderable(renderItem, values, state);
}