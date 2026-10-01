/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type Renderable, type RenderableState, createRenderable } from '../renderable.js';
import type { WebGLContext } from '../webgl/context.js';
import { createGraphicsRenderItem, type Transparency } from '../webgl/render-item.js';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, DefineSpec, type Values, InternalSchema, SizeSchema, type InternalValues, GlobalTextureSchema, type GlobalDefineValues, type GlobalDefines, GlobalDefineSchema, AnimationSchema } from './schema.js';
import { PointsShaderCode } from '../shader-code.js';
import { ValueCell } from '@molstar/core/util';

export const PointsSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aGroup: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),
    dPointSizeAttenuation: DefineSpec('boolean'),
    dPointStyle: DefineSpec('string', ['square', 'circle', 'fuzzy']),

    ...AnimationSchema,
};
export type PointsSchema = typeof PointsSchema
export type PointsValues = Values<PointsSchema>

export function PointsRenderable(ctx: WebGLContext, id: number, values: PointsValues, state: RenderableState, materialId: number, transparency: Transparency, globals: GlobalDefines): Renderable<PointsValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...GlobalDefineSchema, ...InternalSchema, ...PointsSchema };
    const renderValues: PointsValues & InternalValues & GlobalDefineValues = {
        ...values,
        uObjectId: ValueCell.create(id),
        dLightCount: ValueCell.create(globals.dLightCount),
        dColorMarker: ValueCell.create(globals.dColorMarker),
    };
    const shaderCode = PointsShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'points', shaderCode, schema, renderValues, materialId, transparency);
    return createRenderable(renderItem, renderValues, state);
}