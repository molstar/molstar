/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Renderable, RenderableState, createRenderable } from '../renderable';
import { WebGLContext } from '../webgl/context';
import { createGraphicsRenderItem } from '../webgl/render-item';
import { AttributeSpec, Values, UniformSpec, GlobalUniformSchema, InternalSchema, TextureSpec, ValueSpec, ElementsSpec, DefineSpec, InternalValues, GlobalTextureSchema } from './schema';
import { DirectVolumeShaderCode } from '../shader-code';
import { ValueCell } from '../../mol-util';

export const DirectVolumeSchema = {
    uColor: UniformSpec('v3'),
    uColorTexDim: UniformSpec('v2'),
    tColor: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    dColorType: DefineSpec('string', ['uniform', 'attribute', 'instance', 'group', 'groupInstance', 'vertex', 'vertexInstance']),

    uMarkerTexDim: UniformSpec('v2'),
    tMarker: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),

    uOverpaintTexDim: UniformSpec('v2'),
    tOverpaint: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest'),
    dOverpaint: DefineSpec('boolean'),

    uTransparencyTexDim: UniformSpec('v2'),
    tTransparency: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dTransparency: DefineSpec('boolean'),
    dTransparencyVariant: DefineSpec('string', ['single', 'multi']),
    transparencyAverage: ValueSpec('number'),

    dClipObjectCount: DefineSpec('number'),
    dClipVariant: DefineSpec('string', ['instance', 'pixel']),
    uClippingTexDim: UniformSpec('v2'),
    tClipping: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dClipping: DefineSpec('boolean'),

    uVertexCount: UniformSpec('i'),
    uInstanceCount: UniformSpec('i'),
    uGroupCount: UniformSpec('i'),
    uInvariantBoundingSphere: UniformSpec('v4'),

    aInstance: AttributeSpec('float32', 1, 1),
    aTransform: AttributeSpec('float32', 16, 1),

    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    alpha: ValueSpec('number'),

    matrix: ValueSpec('m4'),
    transform: ValueSpec('float32'),
    extraTransform: ValueSpec('float32'),
    hasReflection: ValueSpec('boolean'),

    boundingSphere: ValueSpec('sphere'),
    invariantBoundingSphere: ValueSpec('sphere'),

    aPosition: AttributeSpec('float32', 3, 0),
    elements: ElementsSpec('uint32'),

    uAlpha: UniformSpec('f'),

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