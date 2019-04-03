/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util'
import { Sphere3D } from 'mol-math/geometry'
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { TransformData } from '../transform-data';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { GeometryUtils } from '../geometry';
import { Theme } from 'mol-theme/theme';
import { Color } from 'mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { IsosurfaceValues } from 'mol-gl/renderable/isosurface';
import { calculateTransformBoundingSphere } from 'mol-gl/renderable/util';
import { Texture } from 'mol-gl/webgl/texture';
import { Vec2 } from 'mol-math/linear-algebra';
import { fillSerial } from 'mol-util/array';

export interface Isosurface {
    readonly kind: 'isosurface',

    /** Number of vertices in the isosurface */
    readonly vertexCount: ValueCell<number>,
    /** Number of groups in the isosurface */
    readonly groupCount: ValueCell<number>,

    readonly vertexTexture: ValueCell<Texture>,
    readonly vertexTextureDim: ValueCell<Vec2>,

    /** Normal buffer as array of xyz values for each vertex wrapped in a value cell */
    readonly normalBuffer: ValueCell<Float32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,

    readonly boundingSphere: ValueCell<Sphere3D>,
}

export namespace Isosurface {
    export function create(vertexCount: number, groupCount: number, vertexTexture: Texture, normalBuffer: Float32Array, groupBuffer: Float32Array, boundingSphere: Sphere3D, isosurface?: Isosurface): Isosurface {
        const { width, height } = vertexTexture
        if (isosurface) {
            ValueCell.update(isosurface.vertexCount, vertexCount)
            ValueCell.update(isosurface.groupCount, groupCount)
            ValueCell.update(isosurface.vertexTexture, vertexTexture)
            ValueCell.update(isosurface.vertexTextureDim, Vec2.set(isosurface.vertexTextureDim.ref.value, width, height))
            ValueCell.update(isosurface.normalBuffer, normalBuffer)
            ValueCell.update(isosurface.groupBuffer, groupBuffer)
            ValueCell.update(isosurface.boundingSphere, boundingSphere)
            return isosurface
        } else {
            return {
                kind: 'isosurface',
                vertexCount: ValueCell.create(vertexCount),
                groupCount: ValueCell.create(groupCount),
                vertexTexture: ValueCell.create(vertexTexture),
                vertexTextureDim: ValueCell.create(Vec2.create(width, height)),
                normalBuffer: ValueCell.create(normalBuffer),
                groupBuffer: ValueCell.create(groupBuffer),
                boundingSphere: ValueCell.create(boundingSphere),
            }
        }
    }

    export function createEmpty(isosurface?: Isosurface): Isosurface {
        return {} as Isosurface // TODO
    }

    export const Params = {
        ...BaseGeometry.Params,
        doubleSided: PD.Boolean(false),
        flipSided: PD.Boolean(false),
        flatShaded: PD.Boolean(false),
    }
    export type Params = typeof Params

    export const Utils: GeometryUtils<Isosurface, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState: BaseGeometry.createRenderableState,
        updateRenderableState: BaseGeometry.updateRenderableState
    }

    function createValues(isosurface: Isosurface, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): IsosurfaceValues {
        const { instanceCount, groupCount } = locationIt
        const color = createColors(locationIt, theme.color)
        const marker = createMarkers(instanceCount * groupCount)
        const overpaint = createEmptyOverpaint()
        const transparency = createEmptyTransparency()

        const counts = { drawCount: isosurface.vertexCount.ref.value, groupCount, instanceCount }

        const transformBoundingSphere = calculateTransformBoundingSphere(isosurface.boundingSphere.ref.value, transform.aTransform.ref.value, transform.instanceCount.ref.value)

        return {
            tPosition: isosurface.vertexTexture,
            uPositionTexDim: isosurface.vertexTextureDim,
            aIndex: ValueCell.create(fillSerial(new Float32Array(isosurface.vertexCount.ref.value))),
            aNormal: isosurface.normalBuffer,
            aGroup: isosurface.groupBuffer,
            boundingSphere: ValueCell.create(transformBoundingSphere),
            invariantBoundingSphere: isosurface.boundingSphere,

            ...color,
            ...marker,
            ...overpaint,
            ...transparency,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            dDoubleSided: ValueCell.create(props.doubleSided),
            dFlatShaded: ValueCell.create(props.flatShaded),
            dFlipSided: ValueCell.create(props.flipSided),
            dPositionTexture: ValueCell.create(true),
        }
    }

    function createValuesSimple(isosurface: Isosurface, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform)
        const p = { ...PD.getDefaultValues(Params), ...props }
        return createValues(isosurface, s.transform, s.locationIterator, s.theme, p)
    }

    function updateValues(values: IsosurfaceValues, props: PD.Values<Params>) {
        if (Color.fromNormalizedArray(values.uHighlightColor.ref.value, 0) !== props.highlightColor) {
            ValueCell.update(values.uHighlightColor, Color.toArrayNormalized(props.highlightColor, values.uHighlightColor.ref.value, 0))
        }
        if (Color.fromNormalizedArray(values.uSelectColor.ref.value, 0) !== props.selectColor) {
            ValueCell.update(values.uSelectColor, Color.toArrayNormalized(props.selectColor, values.uSelectColor.ref.value, 0))
        }
        ValueCell.updateIfChanged(values.alpha, props.alpha) // `uAlpha` is set in renderable.render
        ValueCell.updateIfChanged(values.dUseFog, props.useFog)

        ValueCell.updateIfChanged(values.dDoubleSided, props.doubleSided)
        ValueCell.updateIfChanged(values.dFlatShaded, props.flatShaded)
        ValueCell.updateIfChanged(values.dFlipSided, props.flipSided)
    }

    function updateBoundingSphere(values: IsosurfaceValues, isosurface: Isosurface) {
        const invariantBoundingSphere = isosurface.boundingSphere.ref.value
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value)
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere)
        }
    }
}