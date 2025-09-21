/**
 * Copyright (c) 2019-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem, Transparency } from '../webgl/render-item';
import { GlobalUniformSchema, BaseSchema, AttributeSpec, UniformSpec, Values, InternalSchema, SizeSchema, InternalValues, TextureSpec, ElementsSpec, ValueSpec, GlobalTextureSchema, GlobalDefineValues, GlobalDefines, GlobalDefineSchema } from './schema';
import { TextShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const TextSchema = {
    ...BaseSchema,
    ...SizeSchema,
    aGroup: AttributeSpec('float32', 1, 0),
    aPosition: AttributeSpec('float32', 3, 0),
    aMapping: AttributeSpec('float32', 2, 0),
    aDepth: AttributeSpec('float32', 1, 0),
    elements: ElementsSpec('uint32'),

    aTexCoord: AttributeSpec('float32', 2, 0),
    tFont: TextureSpec('image-uint8', 'alpha', 'ubyte', 'linear'),
    padding: ValueSpec('number'),

    uBorderWidth: UniformSpec('f', 'material'),
    uBorderColor: UniformSpec('v3', 'material'),
    uOffsetX: UniformSpec('f', 'material'),
    uOffsetY: UniformSpec('f', 'material'),
    uOffsetZ: UniformSpec('f', 'material'),
    uBackgroundColor: UniformSpec('v3', 'material'),
    uBackgroundOpacity: UniformSpec('f', 'material'),
};
export type TextSchema = typeof TextSchema
export type TextValues = Values<TextSchema>

export function TextRenderable(ctx: WebGLContext, id: number, values: TextValues, state: RenderableState, materialId: number, transparency: Transparency, globals: GlobalDefines): Renderable<TextValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...GlobalDefineSchema, ...InternalSchema, ...TextSchema };
    const renderValues: TextValues & InternalValues & GlobalDefineValues = {
        ...values,
        uObjectId: ValueCell.create(id),
        dLightCount: ValueCell.create(globals.dLightCount),
        dColorMarker: ValueCell.create(globals.dColorMarker),
    };
    const shaderCode = TextShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, renderValues, materialId, transparency);
    return createRenderable(renderItem, renderValues, state);
}