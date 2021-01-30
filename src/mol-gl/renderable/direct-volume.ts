/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, GlobalUniformSchema, InternalSchema, TextureSpec, ElementsSpec, DefineSpec, InternalValues, GlobalTextureSchema, BaseSchema } from './schema';
import { DirectVolumeShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const DirectVolumeSchema = {
    ...BaseSchema,

    aPosition: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),

    uColor: UniformSpec('v3'),
    uColorTexDim: UniformSpec('v2'),
    tColor: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    dColorType: DefineSpec('string', ['uniform', 'attribute', 'instance', 'group', 'groupInstance', 'vertex', 'vertexInstance']),

    uIsoValue: UniformSpec('v2'),
    uBboxMin: UniformSpec('v3'),
    uBboxMax: UniformSpec('v3'),
    uBboxSize: UniformSpec('v3'),
    uMaxSteps: UniformSpec('i'),
    uStepScale: UniformSpec('f'),
    uJumpLength: UniformSpec('f'),
    uTransform: UniformSpec('m4'),
    uGridDim: UniformSpec('v3'),
    dRenderMode: DefineSpec('string', ['isosurface', 'volume']),
    dSingleLayer: DefineSpec('boolean'),
    tTransferTex: TextureSpec('image-uint8', 'rgba', 'ubyte', 'linear'),
    uTransferScale: UniformSpec('f'),

    dGridTexType: DefineSpec('string', ['2d', '3d']),
    uGridTexDim: UniformSpec('v3'),
    tGridTex: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    uGridStats: UniformSpec('v4'), // [min, max, mean, sigma]

    uCellDim: UniformSpec('v3'),
    uCartnToUnit: UniformSpec('m4'),
    uUnitToCartn: UniformSpec('m4'),
    dPackedGroup: DefineSpec('boolean'),

    dDoubleSided: DefineSpec('boolean'),
    dFlipSided: DefineSpec('boolean'),
    dFlatShaded: DefineSpec('boolean'),
    dIgnoreLight: DefineSpec('boolean'),
    dXrayShaded: DefineSpec('boolean'),
};
export type DirectVolumeSchema = typeof DirectVolumeSchema
export type DirectVolumeValues = Values<DirectVolumeSchema>

export function DirectVolumeRenderable(ctx: WebGLContext, id: number, values: DirectVolumeValues, state: RenderableState, materialId: number): Renderable<DirectVolumeValues> {
    const schema = { ...GlobalUniformSchema, ...GlobalTextureSchema, ...InternalSchema, ...DirectVolumeSchema };
    if (!ctx.isWebGL2) {
        // workaround for webgl1 limitation that loop counters need to be `const`
        (schema.uMaxSteps as any) = DefineSpec('number');
    }
    const internalValues: InternalValues = {
        uObjectId: ValueCell.create(id),
    };
    const shaderCode = DirectVolumeShaderCode;
    const renderItem = createGraphicsRenderItem(ctx, 'triangles', shaderCode, schema, { ...values, ...internalValues }, materialId);
    return createRenderable(renderItem, values, state);
}