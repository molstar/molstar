/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util'
import { Mat4 } from 'mol-math/linear-algebra'
import { transformPositionArray/* , transformDirectionArray, getNormalMatrix */ } from '../../util';
import { GeometryUtils } from '../geometry';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { createSizes } from '../size-data';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../util/location-iterator';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { calculateBoundingSphere } from 'mol-gl/renderable/util';
import { Sphere3D } from 'mol-math/geometry';
import { Theme } from 'mol-theme/theme';
import { PointsValues } from 'mol-gl/renderable/points';
import { RenderableState } from 'mol-gl/renderable';
import { Color } from 'mol-util/color';
import { BaseGeometry } from '../base';

/** Point cloud */
export interface Points {
    readonly kind: 'points',
    /** Number of vertices in the point cloud */
    pointCount: number,
    /** Center buffer as array of xyz values wrapped in a value cell */
    readonly centerBuffer: ValueCell<Float32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,
}

export namespace Points {
    export function createEmpty(points?: Points): Points {
        const cb = points ? points.centerBuffer.ref.value : new Float32Array(0)
        const gb = points ? points.groupBuffer.ref.value : new Float32Array(0)
        return {
            kind: 'points',
            pointCount: 0,
            centerBuffer: points ? ValueCell.update(points.centerBuffer, cb) : ValueCell.create(cb),
            groupBuffer: points ? ValueCell.update(points.groupBuffer, gb) : ValueCell.create(gb),
        }
    }

    export function transformImmediate(points: Points, t: Mat4) {
        transformRangeImmediate(points, t, 0, points.pointCount)
    }

    export function transformRangeImmediate(points: Points, t: Mat4, offset: number, count: number) {
        const c = points.centerBuffer.ref.value
        transformPositionArray(t, c, offset, count)
        ValueCell.update(points.centerBuffer, c);
    }

    //

    export const Params = {
        ...BaseGeometry.Params,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        pointSizeAttenuation: PD.Boolean(false),
        pointFilledCircle: PD.Boolean(false),
        pointEdgeBleach: PD.Numeric(0.2, { min: 0, max: 1, step: 0.05 }),
    }
    export type Params = typeof Params

    export const Utils: GeometryUtils<Points, Params> = {
        Params,
        createEmpty,
        createValues,
        createValuesSimple,
        updateValues,
        updateBoundingSphere,
        createRenderableState,
        updateRenderableState
    }

    function createValues(points: Points, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): PointsValues {
        const { instanceCount, groupCount } = locationIt
        const color = createColors(locationIt, theme.color)
        const size = createSizes(locationIt, theme.size)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: points.pointCount, groupCount, instanceCount }

        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            points.centerBuffer.ref.value, points.pointCount,
            transform.aTransform.ref.value, transform.instanceCount.ref.value
        )

        return {
            aPosition: points.centerBuffer,
            aGroup: points.groupBuffer,
            boundingSphere: ValueCell.create(boundingSphere),
            invariantBoundingSphere: ValueCell.create(invariantBoundingSphere),
            ...color,
            ...size,
            ...marker,
            ...transform,

            ...BaseGeometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor),
            dPointSizeAttenuation: ValueCell.create(props.pointSizeAttenuation),
            dPointFilledCircle: ValueCell.create(props.pointFilledCircle),
            uPointEdgeBleach: ValueCell.create(props.pointEdgeBleach),
        }
    }

    function createValuesSimple(points: Points, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform)
        const p = { ...PD.getDefaultValues(Params), ...props }
        return createValues(points, s.transform, s.locationIterator, s.theme, p)
    }

    function updateValues(values: PointsValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props)
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor)
        ValueCell.updateIfChanged(values.dPointSizeAttenuation, props.pointSizeAttenuation)
        ValueCell.updateIfChanged(values.dPointFilledCircle, props.pointFilledCircle)
        ValueCell.updateIfChanged(values.uPointEdgeBleach, props.pointEdgeBleach)
    }

    function updateBoundingSphere(values: PointsValues, points: Points) {
        const { boundingSphere, invariantBoundingSphere } = calculateBoundingSphere(
            values.aPosition.ref.value, points.pointCount,
            values.aTransform.ref.value, values.instanceCount.ref.value
        )
        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere)
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere)
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props)
        updateRenderableState(state, props)
        return state
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props)
        state.opaque = state.opaque && (
            !props.pointFilledCircle ||
            (props.pointFilledCircle && props.pointEdgeBleach === 0)
        )
    }
}