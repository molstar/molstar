/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '../../../mol-util';
import { Mat4, Vec4 } from '../../../mol-math/linear-algebra';
import { transformPositionArray, GroupMapping, createGroupMapping} from '../../util';
import { GeometryUtils } from '../geometry';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { createSizes } from '../size-data';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../util/location-iterator';
import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { calculateInvariantBoundingSphere, calculateTransformBoundingSphere } from '../../../mol-gl/renderable/util';
import { Sphere3D } from '../../../mol-math/geometry';
import { Theme } from '../../../mol-theme/theme';
import { PointsValues } from '../../../mol-gl/renderable/points';
import { RenderableState } from '../../../mol-gl/renderable';
import { Color } from '../../../mol-util/color';
import { BaseGeometry } from '../base';
import { createEmptyOverpaint } from '../overpaint-data';
import { createEmptyTransparency } from '../transparency-data';
import { hashFnv32a } from '../../../mol-data/util';
import { createEmptyClipping } from '../clipping-data';

/** Point cloud */
export interface Points {
    readonly kind: 'points',

    /** Number of vertices in the point cloud */
    pointCount: number,

    /** Center buffer as array of xyz values wrapped in a value cell */
    readonly centerBuffer: ValueCell<Float32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,

    /** Bounding sphere of the points */
    readonly boundingSphere: Sphere3D
    /** Maps group ids to point indices */
    readonly groupMapping: GroupMapping

    setBoundingSphere(boundingSphere: Sphere3D): void
}

export namespace Points {
    export function create(centers: Float32Array, groups: Float32Array, pointCount: number, points?: Points): Points {
        return points ?
            update(centers, groups, pointCount, points) :
            fromArrays(centers, groups, pointCount);
    }

    export function createEmpty(points?: Points): Points {
        const cb = points ? points.centerBuffer.ref.value : new Float32Array(0);
        const gb = points ? points.groupBuffer.ref.value : new Float32Array(0);
        return create(cb, gb, 0, points);
    }

    function hashCode(points: Points) {
        return hashFnv32a([
            points.pointCount, points.centerBuffer.ref.version, points.groupBuffer.ref.version,
        ]);
    }

    function fromArrays(centers: Float32Array, groups: Float32Array, pointCount: number): Points {

        const boundingSphere = Sphere3D();
        let groupMapping: GroupMapping;

        let currentHash = -1;
        let currentGroup = -1;

        const points = {
            kind: 'points' as const,
            pointCount,
            centerBuffer: ValueCell.create(centers),
            groupBuffer: ValueCell.create(groups),
            get boundingSphere() {
                const newHash = hashCode(points);
                if (newHash !== currentHash) {
                    const b = calculateInvariantBoundingSphere(points.centerBuffer.ref.value, points.pointCount, 1);
                    Sphere3D.copy(boundingSphere, b);
                    currentHash = newHash;
                }
                return boundingSphere;
            },
            get groupMapping() {
                if (points.groupBuffer.ref.version !== currentGroup) {
                    groupMapping = createGroupMapping(points.groupBuffer.ref.value, points.pointCount);
                    currentGroup = points.groupBuffer.ref.version;
                }
                return groupMapping;
            },
            setBoundingSphere(sphere: Sphere3D) {
                Sphere3D.copy(boundingSphere, sphere);
                currentHash = hashCode(points);
            }
        };
        return points;
    }

    function update(centers: Float32Array, groups: Float32Array, pointCount: number, points: Points) {
        points.pointCount = pointCount;
        ValueCell.update(points.centerBuffer, centers);
        ValueCell.update(points.groupBuffer, groups);
        return points;
    }

    export function transform(points: Points, t: Mat4) {
        const c = points.centerBuffer.ref.value;
        transformPositionArray(t, c, 0, points.pointCount);
        ValueCell.update(points.centerBuffer, c);
    }

    //

    export const Params = {
        ...BaseGeometry.Params,
        sizeFactor: PD.Numeric(1, { min: 0, max: 10, step: 0.1 }),
        pointSizeAttenuation: PD.Boolean(false),
        pointFilledCircle: PD.Boolean(false),
        pointEdgeBleach: PD.Numeric(0.2, { min: 0, max: 1, step: 0.05 }),
    };
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
    };

    function createValues(points: Points, transform: TransformData, locationIt: LocationIterator, theme: Theme, props: PD.Values<Params>): PointsValues {
        const { instanceCount, groupCount } = locationIt;
        const color = createColors(locationIt, theme.color);
        const size = createSizes(locationIt, theme.size);
        const marker = createMarkers(instanceCount * groupCount);
        const overpaint = createEmptyOverpaint();
        const transparency = createEmptyTransparency();
        const clipping = createEmptyClipping();

        const counts = { drawCount: points.pointCount, groupCount, instanceCount };

        const invariantBoundingSphere = Sphere3D.clone(points.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, transform.aTransform.ref.value, instanceCount);

        return {
            aPosition: points.centerBuffer,
            aGroup: points.groupBuffer,
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

            ...BaseGeometry.createValues(props, counts),
            uSizeFactor: ValueCell.create(props.sizeFactor),
            dPointSizeAttenuation: ValueCell.create(props.pointSizeAttenuation),
            dPointFilledCircle: ValueCell.create(props.pointFilledCircle),
            uPointEdgeBleach: ValueCell.create(props.pointEdgeBleach),
        };
    }

    function createValuesSimple(points: Points, props: Partial<PD.Values<Params>>, colorValue: Color, sizeValue: number, transform?: TransformData) {
        const s = BaseGeometry.createSimple(colorValue, sizeValue, transform);
        const p = { ...PD.getDefaultValues(Params), ...props };
        return createValues(points, s.transform, s.locationIterator, s.theme, p);
    }

    function updateValues(values: PointsValues, props: PD.Values<Params>) {
        BaseGeometry.updateValues(values, props);
        ValueCell.updateIfChanged(values.uSizeFactor, props.sizeFactor);
        ValueCell.updateIfChanged(values.dPointSizeAttenuation, props.pointSizeAttenuation);
        ValueCell.updateIfChanged(values.dPointFilledCircle, props.pointFilledCircle);
        ValueCell.updateIfChanged(values.uPointEdgeBleach, props.pointEdgeBleach);
    }

    function updateBoundingSphere(values: PointsValues, points: Points) {
        const invariantBoundingSphere = Sphere3D.clone(points.boundingSphere);
        const boundingSphere = calculateTransformBoundingSphere(invariantBoundingSphere, values.aTransform.ref.value, values.instanceCount.ref.value);

        if (!Sphere3D.equals(boundingSphere, values.boundingSphere.ref.value)) {
            ValueCell.update(values.boundingSphere, boundingSphere);
        }
        if (!Sphere3D.equals(invariantBoundingSphere, values.invariantBoundingSphere.ref.value)) {
            ValueCell.update(values.invariantBoundingSphere, invariantBoundingSphere);
            ValueCell.update(values.uInvariantBoundingSphere, Vec4.fromSphere(values.uInvariantBoundingSphere.ref.value, invariantBoundingSphere));
        }
    }

    function createRenderableState(props: PD.Values<Params>): RenderableState {
        const state = BaseGeometry.createRenderableState(props);
        updateRenderableState(state, props);
        return state;
    }

    function updateRenderableState(state: RenderableState, props: PD.Values<Params>) {
        BaseGeometry.updateRenderableState(state, props);
        state.opaque = state.opaque && (
            !props.pointFilledCircle ||
            (props.pointFilledCircle && props.pointEdgeBleach === 0)
        );
        state.writeDepth = state.opaque;
    }
}