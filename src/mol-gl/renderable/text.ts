/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, Values, InternalSchema, SizeSchema, InternalValues, TextureSpec, ElementsSpec, ValueSpec } from './schema';
import { TextShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const TextSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aPosition: AttributeSpec('float32', 3, 0),
    aMapping: AttributeSpec('float32', 2, 0),
    aDepth: AttributeSpec('float32', 1, 0),
    elements: ElementsSpec('uint32'),

    aTexCoord: AttributeSpec('float32', 2, 0),
    tFont: TextureSpec('image-uint8', 'alpha', 'ubyte', 'linear'),
    padding: ValueSpec('number'),

    uBorderWidth: UniformSpec('f'),
    uBorderColor: UniformSpec('v3'),
    uOffsetX: UniformSpec('f'),
    uOffsetY: UniformSpec('f'),
    uOffsetZ: UniformSpec('f'),
    uBackgroundColor: UniformSpec('v3'),
    uBackgroundOpacity: UniformSpec('f'),
};
export type TextSchema = typeof TextSchema
export type TextValues = Values<TextSchema>

export function TextRenderable(ctx: WebGLContext, id: number, values: TextValues, state: RenderableState, materialId: number): Renderable<TextValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...TextSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = TextShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId);
    return createRenderable(renderItem, values, state);
}