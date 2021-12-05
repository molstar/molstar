/**
 * Copyright (c) 2018-2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../mol-util';
import { AttributeItemSize, ElementsKind, AttributeValues, AttributeKind, DataTypeArrayType } from '../webgl/buffer';
import { UniformKind, UniformValues, UniformKindValue } from '../webgl/uniform';
import { DefineKind, DefineValues } from '../shader-code';
import { Mat4 } from '../../mol-math/linear-algebra';
import { TextureValues, TextureType, TextureFormat, TextureFilter, TextureKind, TextureKindValue } from '../webgl/texture';
import { Sphere3D } from '../../mol-math/geometry';

export type ValueKindType = {
    'number': number
    'string': string
    'boolean': boolean
    'any': any
    'unknown': unknown

    'm4': Mat4,
    'float32': Float32Array
    'sphere': Sphere3D
}
export type ValueKind = keyof ValueKindType

//

export type KindValue = UniformKindValue & DataTypeArrayType & TextureKindValue & ValueKindType

export type Values<S extends RenderableSchema> = { readonly [k in keyof S]: ValueCell<KindValue[S[k]['kind']]> }
export type UnboxedValues<S extends RenderableSchema> = { readonly [k in keyof S]: KindValue[S[k]['kind']] }

export function splitValues(schema: RenderableSchema, values: RenderableValues) {
    const attributeValues: AttributeValues = {};
    const defineValues: DefineValues = {};
    const textureValues: TextureValues = {};
    const uniformValues: UniformValues = {};
    const materialUniformValues: UniformValues = {};
    const bufferedUniformValues: UniformValues = {};
    Object.keys(schema).forEach(k => {
        const spec = schema[k];
        if (spec.type === 'attribute') attributeValues[k] = values[k];
        if (spec.type === 'define') defineValues[k] = values[k];
        // check if k exists in values to exclude global textures
        if (spec.type === 'texture' && values[k] !== undefined) textureValues[k] = values[k];
        // check if k exists in values to exclude global uniforms
        if (spec.type === 'uniform' && values[k] !== undefined) {
            if (spec.variant === 'material') materialUniformValues[k] = values[k];
            else if (spec.variant === 'buffered') bufferedUniformValues[k] = values[k];
            else uniformValues[k] = values[k];
        }
    });
    return { attributeValues, defineValues, textureValues, uniformValues, materialUniformValues, bufferedUniformValues };
}

export type Versions<T extends RenderableValues> = { [k in keyof T]: number }
export function getValueVersions<T extends RenderableValues>(values: T) {
    const versions: Versions<any> = {};
    Object.keys(values).forEach(k => {
        versions[k] = values[k].ref.version;
    });
    return versions as Versions<T>;
}

//

export type AttributeSpec<K extends AttributeKind> = { type: 'attribute', kind: K, itemSize: AttributeItemSize, divisor: number }
export function AttributeSpec<K extends AttributeKind>(kind: K, itemSize: AttributeItemSize, divisor: number): AttributeSpec<K> {
    return { type: 'attribute', kind, itemSize, divisor };
}

export type UniformSpec<K extends UniformKind> = { type: 'uniform', kind: K, variant?: 'material' | 'buffered' }
export function UniformSpec<K extends UniformKind>(kind: K, variant?: 'material' | 'buffered'): UniformSpec<K> {
    return { type: 'uniform', kind, variant };
}

export type TextureSpec<K extends TextureKind> = { type: 'texture', kind: K, format: TextureFormat, dataType: TextureType, filter: TextureFilter }
export function TextureSpec<K extends TextureKind>(kind: K, format: TextureFormat, dataType: TextureType, filter: TextureFilter): TextureSpec<K> {
    return { type: 'texture', kind, format, dataType, filter };
}

export type ElementsSpec<K extends ElementsKind> = { type: 'elements', kind: K }
export function ElementsSpec<K extends ElementsKind>(kind: K): ElementsSpec<K> {
    return { type: 'elements', kind };
}

export type DefineSpec<K extends DefineKind> = { type: 'define', kind: K, options?: string[] }
export function DefineSpec<K extends DefineKind>(kind: K, options?: string[]): DefineSpec<K> {
    return { type: 'define', kind, options };
}

export type ValueSpec<K extends ValueKind> = { type: 'value', kind: K }
export function ValueSpec<K extends ValueKind>(kind: K): ValueSpec<K> {
    return { type: 'value', kind };
}

//

export type RenderableSchema = {
    readonly [k: string]: (
        AttributeSpec<AttributeKind> | UniformSpec<UniformKind> | TextureSpec<TextureKind> |
        ValueSpec<ValueKind> | DefineSpec<DefineKind> | ElementsSpec<ElementsKind>
    )
}
export type RenderableValues = { readonly [k: string]: ValueCell<any> }

//

export const GlobalUniformSchema = {
    uModel: UniformSpec('m4'),
    uView: UniformSpec('m4'),
    uInvView: UniformSpec('m4'),
    uModelView: UniformSpec('m4'),
    uInvModelView: UniformSpec('m4'),
    uProjection: UniformSpec('m4'),
    uInvProjection: UniformSpec('m4'),
    uModelViewProjection: UniformSpec('m4'),
    uInvModelViewProjection: UniformSpec('m4'),

    uIsOrtho: UniformSpec('f'),
    uPixelRatio: UniformSpec('f'),
    uViewport: UniformSpec('v4'),
    uViewOffset: UniformSpec('v2'),
    uDrawingBufferSize: UniformSpec('v2'),

    uCameraPosition: UniformSpec('v3'),
    uCameraDir: UniformSpec('v3'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
    uFogColor: UniformSpec('v3'),

    uTransparentBackground: UniformSpec('b'),

    uClipObjectType: UniformSpec('i[]'),
    uClipObjectInvert: UniformSpec('b[]'),
    uClipObjectPosition: UniformSpec('v3[]'),
    uClipObjectRotation: UniformSpec('v4[]'),
    uClipObjectScale: UniformSpec('v3[]'),

    uLightDirection: UniformSpec('v3[]'),
    uLightColor: UniformSpec('v3[]'),
    uAmbientColor: UniformSpec('v3'),

    uPickingAlphaThreshold: UniformSpec('f'),

    uInteriorDarkening: UniformSpec('f'),
    uInteriorColorFlag: UniformSpec('b'),
    uInteriorColor: UniformSpec('v3'),

    uHighlightColor: UniformSpec('v3'),
    uSelectColor: UniformSpec('v3'),
    uHighlightStrength: UniformSpec('f'),
    uSelectStrength: UniformSpec('f'),
    uMarkerPriority: UniformSpec('i'),

    uXrayEdgeFalloff: UniformSpec('f'),

    uRenderWboit: UniformSpec('b'),
    uMarkingDepthTest: UniformSpec('b'),
} as const;
export type GlobalUniformSchema = typeof GlobalUniformSchema
export type GlobalUniformValues = Values<GlobalUniformSchema>

export const GlobalTextureSchema = {
    tDepth: TextureSpec('texture', 'depth', 'ushort', 'nearest'),
} as const;
export type GlobalTextureSchema = typeof GlobalTextureSchema
export type GlobalTextureValues = Values<GlobalTextureSchema>

export const InternalSchema = {
    uObjectId: UniformSpec('i'),
} as const;
export type InternalSchema = typeof InternalSchema
export type InternalValues = Values<InternalSchema>

export const ColorSchema = {
    // aColor: AttributeSpec('float32', 3, 0), // TODO
    uColor: UniformSpec('v3', 'material'),
    uColorTexDim: UniformSpec('v2'),
    uColorGridDim: UniformSpec('v3'),
    uColorGridTransform: UniformSpec('v4'),
    tColor: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    tPalette: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    tColorGrid: TextureSpec('texture', 'rgb', 'ubyte', 'linear'),
    dColorType: DefineSpec('string', ['uniform', 'attribute', 'instance', 'group', 'groupInstance', 'vertex', 'vertexInstance', 'volume', 'volumeInstance']),
    dUsePalette: DefineSpec('boolean'),
} as const;
export type ColorSchema = typeof ColorSchema
export type ColorValues = Values<ColorSchema>

export const SizeSchema = {
    // aSize: AttributeSpec('float32', 1, 0), // TODO
    uSize: UniformSpec('f', 'material'),
    uSizeTexDim: UniformSpec('v2'),
    tSize: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    dSizeType: DefineSpec('string', ['uniform', 'attribute', 'instance', 'group', 'groupInstance']),
    uSizeFactor: UniformSpec('f'),
} as const;
export type SizeSchema = typeof SizeSchema
export type SizeValues = Values<SizeSchema>

export const MarkerSchema = {
    uMarker: UniformSpec('f'),
    uMarkerTexDim: UniformSpec('v2'),
    tMarker: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dMarkerType: DefineSpec('string', ['uniform', 'groupInstance']),
    markerAverage: ValueSpec('number'),
    markerStatus: ValueSpec('number'),
} as const;
export type MarkerSchema = typeof MarkerSchema
export type MarkerValues = Values<MarkerSchema>

export const OverpaintSchema = {
    uOverpaintTexDim: UniformSpec('v2'),
    tOverpaint: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest'),
    dOverpaint: DefineSpec('boolean'),

    uOverpaintGridDim: UniformSpec('v3'),
    uOverpaintGridTransform: UniformSpec('v4'),
    tOverpaintGrid: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    dOverpaintType: DefineSpec('string', ['groupInstance', 'volumeInstance']),
} as const;
export type OverpaintSchema = typeof OverpaintSchema
export type OverpaintValues = Values<OverpaintSchema>

export const TransparencySchema = {
    uTransparencyTexDim: UniformSpec('v2'),
    tTransparency: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dTransparency: DefineSpec('boolean'),
    transparencyAverage: ValueSpec('number'),

    uTransparencyGridDim: UniformSpec('v3'),
    uTransparencyGridTransform: UniformSpec('v4'),
    tTransparencyGrid: TextureSpec('texture', 'alpha', 'ubyte', 'linear'),
    dTransparencyType: DefineSpec('string', ['groupInstance', 'volumeInstance']),
} as const;
export type TransparencySchema = typeof TransparencySchema
export type TransparencyValues = Values<TransparencySchema>

export const SubstanceSchema = {
    uSubstanceTexDim: UniformSpec('v2'),
    tSubstance: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest'),
    dSubstance: DefineSpec('boolean'),

    uSubstanceGridDim: UniformSpec('v3'),
    uSubstanceGridTransform: UniformSpec('v4'),
    tSubstanceGrid: TextureSpec('texture', 'rgba', 'ubyte', 'linear'),
    dSubstanceType: DefineSpec('string', ['groupInstance', 'volumeInstance']),
} as const;
export type SubstanceSchema = typeof SubstanceSchema
export type SubstanceValues = Values<SubstanceSchema>

export const ClippingSchema = {
    dClipObjectCount: DefineSpec('number'),
    dClipVariant: DefineSpec('string', ['instance', 'pixel']),

    uClippingTexDim: UniformSpec('v2'),
    tClipping: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dClipping: DefineSpec('boolean'),
} as const;
export type ClippingSchema = typeof ClippingSchema
export type ClippingValues = Values<ClippingSchema>

export const BaseSchema = {
    ...ColorSchema,
    ...MarkerSchema,
    ...OverpaintSchema,
    ...TransparencySchema,
    ...SubstanceSchema,
    ...ClippingSchema,

    dLightCount: DefineSpec('number'),

    aInstance: AttributeSpec('float32', 1, 1),
    /**
     * final per-instance transform calculated for instance `i` as
     * `aTransform[i] = matrix * transform[i] * extraTransform[i]`
     */
    aTransform: AttributeSpec('float32', 16, 1),

    /**
     * final alpha, calculated as `values.alpha * state.alpha`
     */
    uAlpha: UniformSpec('f', 'material'),
    uMetalness: UniformSpec('f', 'material'),
    uRoughness: UniformSpec('f', 'material'),
    uBumpiness: UniformSpec('f', 'material'),

    uVertexCount: UniformSpec('i'),
    uInstanceCount: UniformSpec('i'),
    uGroupCount: UniformSpec('i'),
    uInvariantBoundingSphere: UniformSpec('v4'),

    drawCount: ValueSpec('number'),
    instanceCount: ValueSpec('number'),

    /** base alpha, see uAlpha  */
    alpha: ValueSpec('number'),

    /** global transform, see aTransform */
    matrix: ValueSpec('m4'),
    /** base per-instance transform, see aTransform */
    transform: ValueSpec('float32'),
    /** additional per-instance transform, see aTransform */
    extraTransform: ValueSpec('float32'),
    /** denotes reflection in transform */
    hasReflection: ValueSpec('boolean'),

    /** bounding sphere taking aTransform into account and encompases all instances */
    boundingSphere: ValueSpec('sphere'),
    /** bounding sphere NOT taking aTransform into account */
    invariantBoundingSphere: ValueSpec('sphere'),
} as const;
export type BaseSchema = typeof BaseSchema
export type BaseValues = Values<BaseSchema>