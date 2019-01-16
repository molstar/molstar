/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from 'mol-util/param-definition';
import { ValueCell } from 'mol-util';
import { Geometry } from '../geometry';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { TransformData, createIdentityTransform } from '../transform-data';
import { Theme } from 'mol-theme/theme';
import { createColors } from '../color-data';
import { createSizes, getMaxSize } from '../size-data';
import { createMarkers } from '../marker-data';
import { ColorNames } from 'mol-util/color/tables';
import { NullLocation } from 'mol-model/location';
import { UniformColorTheme } from 'mol-theme/color/uniform';
import { UniformSizeTheme } from 'mol-theme/size/uniform';
import { Sphere3D } from 'mol-math/geometry';
import { calculateBoundingSphere } from 'mol-gl/renderable/util';
import { TextValues } from 'mol-gl/renderable/text';
import { Color } from 'mol-util/color';
import { Vec3 } from 'mol-math/linear-algebra';
import { FontAtlas, getFontAtlas, FontAtlasParams } from './font-atlas';
import { RenderableState } from 'mol-gl/renderable';
import { clamp } from 'mol-math/interpolate';

type TextAttachment = 'bottom-left' | 'bottom-center' | 'bottom-right' | 'middle-left' | 'middle-center' | 'middle-right' | 'top-left' | 'top-center' | 'top-right'

/** Text */
export interface Text {
    readonly kind: 'text',

    /** Number of characters in the text */
    readonly charCount: number,
    /** Font Atlas */
    readonly fontAtlas: FontAtlas,

    /** Center buffer as array of xyz values wrapped in a value cell */
    readonly centerBuffer: ValueCell<Float32Array>,
    /** Mapping buffer as array of xy values wrapped in a value cell */
    readonly mappingBuffer: ValueCell<Float32Array>,
    /** Index buffer as array of center index triplets wrapped in a value cell */
    readonly indexBuffer: ValueCell<Uint32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,
    /** Texture coordinates buffer as array of uv values wrapped in a value cell */
    readonly tcoordBuffer: ValueCell<Float32Array>,
}

export namespace Text {
    export function createEmpty(text?: Text): Text {
        const cb = text ? text.centerBuffer.ref.value : new Float32Array(0)
        const mb = text ? text.mappingBuffer.ref.value : new Float32Array(0)
        const ib = text ? text.indexBuffer.ref.value : new Uint32Array(0)
        const gb = text ? text.groupBuffer.ref.value : new Float32Array(0)
        const tb = text ? text.tcoordBuffer.ref.value : new Float32Array(0)
        return {
            kind: 'text',
            charCount: 0,
            fontAtlas: getFontAtlas({}),
            centerBuffer: text ? ValueCell.update(text.centerBuffer, cb) : ValueCell.create(cb),
            mappingBuffer: text ? ValueCell.update(text.mappingBuffer, mb) : ValueCell.create(mb),
            indexBuffer: text ? ValueCell.update(text.indexBuffer, ib) : ValueCell.create(ib),
            groupBuffer: text ? ValueCell.update(text.groupBuffer, gb) : ValueCell.create(gb),
            tcoordBuffer: text ? ValueCell.update(text.tcoordBuffer, tb) : ValueCell.create(tb)
        }
    }

    export const Params = {
        ...Geometry.Params,
        ...FontAtlasParams,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),

        borderWidth: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),
        borderColor: PD.Color(ColorNames.grey),
        offsetX: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        offsetY: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        offsetZ: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        background: PD.Boolean(false),
        backgroundMargin: PD.Numeric(0.2, { min: 0, max: 10, step: 0.1 }),
        backgroundColor: PD.Color(ColorNames.grey),
        backgroundOpacity: PD.Numeric(0, { min: 0, max: 1, step: 0.01 }),

        attachment: PD.Select('normal', [['bottom-left', 'bottom-left'], ['bottom-center', 'bottom-center'], ['bottom-right', 'bottom-right'], ['middle-left', 'middle-left'], ['top-left', 'top-left'], ['top-center', 'top-center'], ['top-right', 'top-right']] as [TextAttachment, string][]),
    }
    export type Params = typeof Params

    export function createValues(text: Text, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): TextValues {
        const { instanceCount, groupCount } = locationIt
        if (instanceCount !== transform.instanceCount.ref.value) {
            throw new Error('instanceCount values in TransformData and LocationIterator differ')
        }

        const color = createColors(locationIt, theme.color)
        const size = createSizes(locationIt, theme.size)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: text.charCount * 2 * 3, groupCount, instanceCount }

        const padding = getMaxSize(size)
        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            text.centerBuffer.ref.value, text.charCount * 4,
            transform.aTransform.ref.value, instanceCount, padding
        )

        console.log(props.sizeFactor, text.fontAtlas.lineHeight, props.fontSize)

        return {
            aPosition: text.centerBuffer,
            aMapping: text.mappingBuffer,
            aGroup: text.groupBuffer,
            elements: text.indexBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            ...color,
            ...size,
            ...marker,
            ...transform,

            aTexCoord: text.tcoordBuffer,
            tFont: ValueCell.create(text.fontAtlas.texture),
            padding: ValueCell.create(padding),

            ...Geometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor / text.fontAtlas.lineHeight),

            uBorderWidth: ValueCell.create(clamp(props.borderWidth / 2, 0, 0.5)),
            uBorderColor: ValueCell.create(Color.toArrayNormalized(props.borderColor, Vec3.zero(), 0)),
            uOffsetX: ValueCell.create(props.offsetX),
            uOffsetY: ValueCell.create(props.offsetY),
            uOffsetZ: ValueCell.create(props.offsetZ),
            uBackgroundColor: ValueCell.create(Color.toArrayNormalized(props.backgroundColor, Vec3.zero(), 0)),
            uBackgroundOpacity: ValueCell.create(props.backgroundOpacity),
        }
    }

    export function createValuesSimple(text: Text, props: Partial<PD.Values<Params>>, colorValue = ColorNames.grey, sizeValue = 1, transform?: TransformData): TextValues {

        if (!transform) transform = createIdentityTransform()
        const locationIterator = LocationIterator(1, transform.instanceCount.ref.value, () => NullLocation, false, () => false)
        const theme: Theme = {
            color: UniformColorTheme({}, { value: colorValue}),
            size: UniformSizeTheme({}, { value: sizeValue})
        }
        const p = { ...PD.getDefaultValues(Params), ...props }

        return createValues(text, transform, locationIterator, theme, p)
    }

    export function updateValues(values: TextValues, props: PD.Values<Params>) {
        Geometry.updateValues(values, props)
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor)
    }

    export function updateBoundingSphere(values: TextValues, text: Text) {
        const padding = getMaxSize(values)
        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            values.aPosition.ref.value, text.charCount * 4,
            values.aTransform.ref.value, values.instanceCount.ref.value, padding
        )
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere)
        }
    }

    export function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = Geometry.createRenderableState(props)
        updateRenderableState(state, props)
        return state
    }

    export function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        Geometry.updateRenderableState(state, props)
        state.opaque = false
    }
}