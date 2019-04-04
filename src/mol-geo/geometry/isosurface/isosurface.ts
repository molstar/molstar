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

    readonly geoTextureDim: ValueCell<Vec2>,
    readonly vertexTexture: ValueCell<Texture>,
    readonly normalTexture: ValueCell<Texture>,
    readonly groupTexture: ValueCell<Texture>,

    readonly boundingSphere: ValueCell<Sphere3D>,
}

export namespace Isosurface {
    export function create(vertexCount: number, groupCount: number, vertexTexture: Texture, normalTexture: Texture, groupTexture: Texture, boundingSphere: Sphere3D, isosurface?: Isosurface): Isosurface {
        const { width, height } = vertexTexture
        if (isosurface) {
            ValueCell.update(isosurface.vertexCount, vertexCount)
            ValueCell.update(isosurface.groupCount, groupCount)
            ValueCell.update(isosurface.geoTextureDim, Vec2.set(isosurface.geoTextureDim.ref.value, width, height))
            ValueCell.update(isosurface.vertexTexture, vertexTexture)
            ValueCell.update(isosurface.normalTexture, normalTexture)
            ValueCell.update(isosurface.groupTexture, groupTexture)
            ValueCell.update(isosurface.boundingSphere, boundingSphere)
            return isosurface
        } else {
            return {
                kind: 'isosurface',
                vertexCount: ValueCell.create(vertexCount),
                groupCount: ValueCell.create(groupCount),
                geoTextureDim: ValueCell.create(Vec2.create(width, height)),
                vertexTexture: ValueCell.create(vertexTexture),
                normalTexture: ValueCell.create(normalTexture),
                groupTexture: ValueCell.create(groupTexture),
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
            uGeoTexDim: isosurface.geoTextureDim,
            tPosition: isosurface.vertexTexture,
            tNormal: isosurface.normalTexture,
            tGroup: isosurface.groupTexture,

            // aGroup is used as a triangle index here and the group id is retirieved from the tGroup texture
            aGroup: ValueCell.create(fillSerial(new Float32Array(isosurface.vertexCount.ref.value))),
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
            dGeoTexture: ValueCell.create(true),
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