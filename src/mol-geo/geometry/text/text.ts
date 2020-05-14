/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { ValueCell } from '../../../mol-util';
import { GeometryUtils } from '../geometry';
import { LocationIterator } from '../../../mol-geo/util/location-iterator';
import { TransformData } from '../transform-data';
import { Theme } from '../../../mol-theme/theme';
import { createColors } from '../color-data';
import { createSizes, getMaxSize } from '../size-data';
import { createMarkers } from '../marker-data';
import { ColorNames } from '../../../mol-util/color/names';
import { Sphere3D } from '../../../mol-math/geometry';
import { TextureImage, createTextureImage, calculateInvariantBoundingSphere, calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { TextValues } from '../../../mol-gl/renderable/text';
import { Color } from '../../../mol-util/color';
import { Vec3, Vec4 } from '../../../mol-math/linear-algebra';
import { FontAtlasParams } from './font-atlas';
import { RenderableState } from '../../../mol-gl/renderable';
import { clamp } from '../../../mol-math/interpolate';
import { createRenderObject as _createRenderObject } from '../../../mol-gl/render-object';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { hashFnv32a } from '../../../mol-data/util';
import { GroupMapping, createGroupMapping } from '../../util';
import { createEmptyClipping } from '../clipping-data';

type TextAttachment = (
    'bottom-left' | 'bottom-center' | 'bottom-right' |
    'middle-left' | 'middle-center' | 'middle-right' |
    'top-left' | 'top-center' | 'top-right'
)

/** Text */
export interface Text {
    readonly kind: 'text',

    /** Number of characters in the text */
    charCount: number,

    /** Font Atlas */
    readonly fontTexture: ValueCell<TextureImage<Uint8Array>>,

    /** Center buffer as array of xyz values wrapped in a value cell */
    readonly centerBuffer: ValueCell<Float32Array>,
    /** Mapping buffer as array of xy values wrapped in a value cell */
    readonly mappingBuffer: ValueCell<Float32Array>,
    /** Depth buffer as array of z values wrapped in a value cell */
    readonly depthBuffer: ValueCell<Float32Array>,
    /** Index buffer as array of center index triplets wrapped in a value cell */
    readonly indexBuffer: ValueCell<Uint32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,
    /** Texture coordinates buffer as array of uv values wrapped in a value cell */
    readonly tcoordBuffer: ValueCell<Float32Array>,

    /** Bounding sphere of the text */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to text indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void
}

export namespace Text {
    export function create(fontTexture: TextureImage<Uint8Array>, centers: Float32Array, mappings: Float32Array, depths: Float32Array, indices: Uint32Array, groups: Float32Array, tcoords: Float32Array, charCount: number, text?: Text): Text {
        return text ?
            update(fontTexture, centers, mappings, depths, indices, groups, tcoords, charCount, text) :
            fromData(fontTexture, centers, mappings, depths, indices, groups, tcoords, charCount);
    }

    export function createEmpty(text?: Text): Text {
        const ft = text ? text.fontTexture.ref.value : createTextureImage(0, 1, Uint8Array);
        const cb = text ? text.centerBuffer.ref.value : new Float32Array(0);
        const mb = text ? text.mappingBuffer.ref.value : new Float32Array(0);
        const db = text ? text.depthBuffer.ref.value : new Float32Array(0);
        const ib = text ? text.indexBuffer.ref.value : new Uint32Array(0);
        const gb = text ? text.groupBuffer.ref.value : new Float32Array(0);
        const tb = text ? text.tcoordBuffer.ref.value : new Float32Array(0);
        return create(ft, cb, mb, db, ib, gb, tb, 0, text);
    }

    function hashCode(text: Text) {
        return hashFnv32a([
            text.charCount, text.fontTexture.ref.version,
            text.centerBuffer.ref.version, text.mappingBuffer.ref.version,
            text.depthBuffer.ref.version, text.indexBuffer.ref.version,
            text.groupBuffer.ref.version, text.tcoordBuffer.ref.version
        ]);
    }

    function fromData(fontTexture: TextureImage<Uint8Array>, centers: Float32Array, mappings: Float32Array, depths: Float32Array, indices: Uint32Array, groups: Float32Array, tcoords: Float32Array, charCount: number): Text {

        const boundingSphere = Sphere3D();
        let groupMapping: GroupMapping;

        let currentHash = -1;
        let currentGroup = -1;

        const text = {
            kind: 'text' as const,
            charCount,
            fontTexture: ValueCell.create(fontTexture),
            centerBuffer: ValueCell.create(centers),
            mappingBuffer: ValueCell.create(mappings),
            depthBuffer: ValueCell.create(depths),
            indexBuffer: ValueCell.create(indices),
            groupBuffer: ValueCell.create(groups),
            tcoordBuffer: ValueCell.create(tcoords),
            get boundingSphere() {
                const newHash = hashCode(text);
                if (newHash !== currentHash) {
                    const b = calculateInvariantBoundingSphere(text.centerBuffer.ref.value, text.charCount * 4, 4);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            get groupMapping() {
                if (text.groupBuffer.ref.version !== currentGroup) {
                    groupMapping = createGroupMapping(text.groupBuffer.ref.value, text.charCount, 4);
                    currentGroup = text.groupBuffer.ref.version;
                }
                return groupMapping;
            },
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(text);
            }
        };
        return text;
    }

    function update(fontTexture: TextureImage<Uint8Array>, centers: Float32Array, mappings: Float32Array, depths: Float32Array, indices: Uint32Array, groups: Float32Array, tcoords: Float32Array, charCount: number, text: Text) {
        text.charCount = charCount;
        ValueCell.update(text.fontTexture, fontTexture);
        ValueCell.update(text.centerBuffer, centers);
        ValueCell.update(text.mappingBuffer, mappings);
        ValueCell.update(text.depthBuffer, depths);
        ValueCell.update(text.indexBuffer, indices);
        ValueCell.update(text.groupBuffer, groups);
        ValueCell.update(text.tcoordBuffer, tcoords);
        return text;
    }

    export const Params = {
        ...BaseGeometry.Params,
        ...FontAtlasParams,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),

        borderWidth: PD.Numeric(0, { min: 0, max: 0.5, step: 0.01 }),
        borderColor: PD.Color(ColorNames.grey),
        offsetX: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
        offsetY: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
        offsetZ: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
        background: PD.Boolean(false),
        backgroundMargin: PD.Numeric(0.2, { min: 0, max: 1, step: 0.01 }),
        backgroundColor: PD.Color(ColorNames.grey),
        backgroundOpacity: PD.Numeric(1, { min: 0, max: 1, step: 0.01 }),
        tether: PD.Boolean(false),
        tetherLength: PD.Numeric(1, { min: 0, max: 5, step: 0.1 }),
        tetherBaseWidth: PD.Numeric(0.3, { min: 0, max: 1, step: 0.01 }),

        attachment: PD.Select('middle-center', [
            ['bottom-left', 'bottom-left'], ['bottom-center', 'bottom-center'], ['bottom-right', 'bottom-right'],
            ['middle-left', 'middle-left'], ['middle-center', 'middle-center'], ['middle-right', 'middle-right'],
            ['top-left', 'top-left'], ['top-center', 'top-center'], ['top-right', 'top-right'],
        ] as [TextAttachment, string][]),
    };
    export type Params = typeof Params

    export const Utils: GeometryUtils<Text, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState,
    };

    function createValues(text: Text, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): TextValues {
        const { instanceCount, groupCount } = locationIt;
        if (instanceCount !== transform.instanceCount.ref.value) {
            throw new Error('instanceCount values in TransformData and LocationIterator differ');
        }

        const color = createColors(locationIt, theme.color);
        const size = createSizes(locationIt, theme.size);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: text.charCount * 2 * 3, groupCount, instanceCount };

        const padding = getPadding(text.mappingBuffer.ref.value, text.depthBuffer.ref.value, text.charCount, getMaxSize(size));
        const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), text.boundingSphere, padding);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount);

        return {
            aPosition: text.centerBuffer,
            aMapping: text.mappingBuffer,
            aDepth: text.depthBuffer,
            aGroup: text.groupBuffer,
            elements: text.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            uInvariantBoundingSphere: ValueCell.create(Vec4.ofSphere(invariantBoundingSphere)),
            ...color,
            ...size,
            ...marker,
            ...overpaint,
            ...transparency,
            ...clipping,
            ...transform,

            aTexCoord: text.tcoordBuffer,
            tFont: text.fontTexture,
            padding: ValueCell.create(padding),

            ...BaseGeometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor),

            uBorderWidth: ValueCell.create(clamp(props.borderWidth / 2, 0, 0.5)),
            uBorderColor: ValueCell.create(Color.toArrayNormalized(props.borderColor, Vec3.zero(), 0)),
            uOffsetX: ValueCell.create(props.offsetX),
            uOffsetY: ValueCell.create(props.offsetY),
            uOffsetZ: ValueCell.create(props.offsetZ),
            uBackgroundColor: ValueCell.create(Color.toArrayNormalized(props.backgroundColor, Vec3.zero(), 0)),
            uBackgroundOpacity: ValueCell.create(props.backgroundOpacity),
        };
    }

    function createValuesSimple(text: Text, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(text, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: TextValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);

        ValueCell.updateIfChanged(values.uBorderWidth, props.borderWidth);
        if (Color.fromNormalizedArray(values.uBorderColor.ref.value, 0) !== props.borderColor) {
            Color.toArrayNormalized(props.borderColor, values.uBorderColor.ref.value, 0);
            ValueCell.update(values.uBorderColor, values.uBorderColor.ref.value);
        }
        ValueCell.updateIfChanged(values.uOffsetX, props.offsetX);
        ValueCell.updateIfChanged(values.uOffsetY, props.offsetY);
        ValueCell.updateIfChanged(values.uOffsetZ, props.offsetZ);
        if (Color.fromNormalizedArray(values.uBackgroundColor.ref.value, 0) !== props.backgroundColor) {
            Color.toArrayNormalized(props.backgroundColor, values.uBackgroundColor.ref.value, 0);
            ValueCell.update(values.uBackgroundColor, values.uBackgroundColor.ref.value);
        }
        ValueCell.updateIfChanged(values.uBackgroundOpacity, props.backgroundOpacity);
    }

    function updateBoundingSphere(values: TextValues, text: Text) {
        const padding = getPadding(values.aMapping.ref.value, values.aDepth.ref.value, text.charCount, getMaxSize(values));
        const invariantBoundingSphere = Sphere3D.expand(Sphere3D(), text.boundingSphere, padding);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
        ValueCell.update(values.padding, padding);
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props);
        updateRenderableState(state, props);
        return state;
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props);
        state.pickable = false;
        state.opaque = false;
        state.writeDepth = true;
    }
}

function getPadding(mappings: Float32Array, depths: Float32Array, charCount: number, maxSize: number) {
    let maxOffset = 0;
    let maxDepth = 0;
    for (let i = 0, il = charCount * 4; i < il; ++i) {
        const i2 = 2 * i;
        const ox = Math.abs(mappings[i2]);
        if (ox > maxOffset) maxOffset = ox;
        const oy = Math.abs(mappings[i2 + 1]);
        if (oy > maxOffset) maxOffset = oy;
        const d = Math.abs(depths[i]);
        if (d > maxDepth) maxDepth = d;
    }
    // console.log(maxDepth + maxSize, maxDepth, maxSize, maxSize + maxSize * maxOffset, depths)
    return Math.max(maxDepth, maxSize + maxSize * maxOffset);
    // return maxSize + maxSize * maxOffset + maxDepth
}