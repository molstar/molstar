/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, GlobalUniformSchema, InternalSchema, TextureSpec, ElementsSpec, DefineSpec, InternalValues, BaseSchema, UniformSpec } from './schema';
import { ImageShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';
import { InterpolationTypeNames } from '../../mol-geo/geometry/image/image';

export const ImageSchema = {
    ...BaseSchema,

    aPosition: AttributeSpec('float32', 3, 0),
    aUv: AttributeSpec('float32', 2, 0),

    elements: ElementsSpec('uint32'),

    uImageTexDim: UniformSpec('v2'),
    tImageTex: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    tGroupTex: TextureSpec('image-float32', 'alpha', 'float', 'nearest'),

    dInterpolation: DefineSpec('string', InterpolationTypeNames),
};
export type ImageSchema = typeof ImageSchema
export type ImageValues = Values<ImageSchema>

export function ImageRenderable(ctx: WebGLContext, id: number, values: ImageValues, state: RenderableState, materialId: number): Renderable<ImageValues> {
    const schema = { ...GlobalUniformSchema, ...InternalSchema, ...ImageSchema };
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = ImageShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId);
    return createRenderable(renderItem, values, state);
}