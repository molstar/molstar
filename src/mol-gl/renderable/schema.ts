/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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

    'm4': Mat4,
    'float32': Float32Array
    'sphere': Sphere3D
}
export type ValueKind = keyof ValueKindType

//

export type KindValue = UniformKindValue & DataTypeArrayType & TextureKindValue & ValueKindType

export type Values<S extends RenderableSchema> = { readonly [k in keyof S]: ValueCell<KindValue[S[k]['kind']]> }

export function splitValues(schema: RenderableSchema, values: RenderableValues) {
    const attributeValues: AttributeValues = {};
    const defineValues: DefineValues = {};
    const textureValues: TextureValues = {};
    const uniformValues: UniformValues = {};
    const materialUniformValues: UniformValues = {};
    Object.keys(schema).forEach(k => {
        const spec = schema[k];
        if (spec.type === 'attribute') attributeValues[k] = values[k];
        if (spec.type === 'define') defineValues[k] = values[k];
        if (spec.type === 'texture') textureValues[k] = values[k];
        // check if k exists in values so that global uniforms are excluded here
        if (spec.type === 'uniform' && values[k] !== undefined) {
            if (spec.isMaterial) materialUniformValues[k] = values[k];
            else uniformValues[k] = values[k];
        }
    });
    return { attributeValues, defineValues, textureValues, uniformValues, materialUniformValues };
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

export type UniformSpec<K extends UniformKind> = { type: 'uniform', kind: K, isMaterial: boolean }
export function UniformSpec<K extends UniformKind>(kind: K, isMaterial = false): UniformSpec<K> {
    return { type: 'uniform', kind, isMaterial };
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
    uViewportHeight: UniformSpec('f'),
    uViewport: UniformSpec('v4'),
    uViewOffset: UniformSpec('v2'),

    uCameraPosition: UniformSpec('v3'),
    uNear: UniformSpec('f'),
    uFar: UniformSpec('f'),
    uFogNear: UniformSpec('f'),
    uFogFar: UniformSpec('f'),
    uFogColor: UniformSpec('v3'),

    uTransparentBackground: UniformSpec('i'),

    uClipObjectType: UniformSpec('i[]'),
    uClipObjectPosition: UniformSpec('v3[]'),
    uClipObjectRotation: UniformSpec('v4[]'),
    uClipObjectScale: UniformSpec('v3[]'),

    // all the following could in principle be per object
    // as a kind of 'material' parameter set
    // would need to test performance implications
    uLightIntensity: UniformSpec('f'),
    uAmbientIntensity: UniformSpec('f'),

    uMetalness: UniformSpec('f'),
    uRoughness: UniformSpec('f'),
    uReflectivity: UniformSpec('f'),

    uPickingAlphaThreshold: UniformSpec('f'),

    uInteriorDarkening: UniformSpec('f'),
    uInteriorColorFlag: UniformSpec('i'),
    uInteriorColor: UniformSpec('v3'),

    uHighlightColor: UniformSpec('v3'),
    uSelectColor: UniformSpec('v3'),
} as const;
export type GlobalUniformSchema = typeof GlobalUniformSchema
export type GlobalUniformValues = Values<GlobalUniformSchema> // { [k in keyof GlobalUniformSchema]: ValueCell<any> }

export const InternalSchema = {
    uObjectId: UniformSpec('i'),
} as const;
export type InternalSchema = typeof InternalSchema
export type InternalValues = { [k in keyof InternalSchema]: ValueCell<any> }

export const ColorSchema = {
    // aColor: AttributeSpec('float32', 3, 0), // TODO
    uColor: UniformSpec('v3', true),
    uColorTexDim: UniformSpec('v2'),
    tColor: TextureSpec('image-uint8', 'rgb', 'ubyte', 'nearest'),
    dColorType: DefineSpec('string', ['uniform', 'attribute', 'instance', 'group', 'group_instance']),
} as const;
export type ColorSchema = typeof ColorSchema
export type ColorValues = Values<ColorSchema>

export const SizeSchema = {
    // aSize: AttributeSpec('float32', 1, 0), // TODO
    uSize: UniformSpec('f', true),
    uSizeTexDim: UniformSpec('v2'),
    tSize: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dSizeType: DefineSpec('string', ['uniform', 'attribute', 'instance', 'group', 'group_instance']),
    uSizeFactor: UniformSpec('f'),
} as const;
export type SizeSchema = typeof SizeSchema
export type SizeValues = Values<SizeSchema>

export const MarkerSchema = {
    uMarkerTexDim: UniformSpec('v2'),
    tMarker: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
} as const;
export type MarkerSchema = typeof MarkerSchema
export type MarkerValues = Values<MarkerSchema>

export const OverpaintSchema = {
    uOverpaintTexDim: UniformSpec('v2'),
    tOverpaint: TextureSpec('image-uint8', 'rgba', 'ubyte', 'nearest'),
    dOverpaint: DefineSpec('boolean'),
} as const;
export type OverpaintSchema = typeof OverpaintSchema
export type OverpaintValues = Values<OverpaintSchema>

export const TransparencySchema = {
    uTransparencyTexDim: UniformSpec('v2'),
    tTransparency: TextureSpec('image-uint8', 'alpha', 'ubyte', 'nearest'),
    dTransparency: DefineSpec('boolean'),
    dTransparencyVariant: DefineSpec('string', ['single', 'multi']),
} as const;
export type TransparencySchema = typeof TransparencySchema
export type TransparencyValues = Values<TransparencySchema>

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
    ...ClippingSchema,

    aInstance: AttributeSpec('float32', 1, 1),
    aGroup: AttributeSpec('float32', 1, 0),
    /**
     * final per-instance transform calculated for instance `i` as
     * `aTransform[i] = matrix * transform[i] * extraTransform[i]`
     */
    aTransform: AttributeSpec('float32', 16, 1),

    /**
     * final alpha, calculated as `values.alpha * state.alpha`
     */
    uAlpha: UniformSpec('f', true),
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

    /** bounding sphere taking aTransform into account and encompases all instances */
    boundingSphere: ValueSpec('sphere'),
    /** bounding sphere NOT taking aTransform into account */
    invariantBoundingSphere: ValueSpec('sphere'),
} as const;
export type BaseSchema = typeof BaseSchema
export type BaseValues = Values<BaseSchema>