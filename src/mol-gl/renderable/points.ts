/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, DefineSpec, Values, InternalSchema, SizeSchema, InternalValues } from './schema';
import { PointsShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const PointsSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aPosition: AttributeSpec('float32', 3, 0),
    dPointSizeAttenuation: DefineSpec('boolean'),
    dPointFilledCircle: DefineSpec('boolean'),
    uPointEdgeBleach: UniformSpec('f'),
};
export type PointsSchema = typeof PointsSchema
export type PointsValues = Values<PointsSchema>

export function PointsRenderable(ctx: WebGLContext, id: number, values: PointsValues, state: RenderableState, materialId: number): Renderable<PointsValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...PointsSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = PointsShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'points', shaderCode, schema, { ...values, ...internalValues }, materialId);
    return createRenderable(renderItem, values, state);
}