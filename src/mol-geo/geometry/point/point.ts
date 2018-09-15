/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util'
import { Mat4 } from 'mol-math/linear-algebra'
import { transformPositionArray/* , transformDirectionArray, getNormalMatrix */ } from '../../util';
import { Geometry } from '../geometry';
import { PointValues } from 'mol-gl/renderable';
import { RuntimeContext } from 'mol-task';
import { createColors } from '../color-data';
import { createMarkers } from '../marker-data';
import { createSizes } from '../size-data';
import { TransformData } from '../transform-data';
import { LocationIterator } from '../../util/location-iterator';
import { SizeThemeProps } from 'mol-view/theme/size';

/** Point cloud */
export interface Point {
    readonly kind: 'point',
    /** Number of vertices in the point cloud */
    vertexCount: number,
    /** Vertex buffer as array of xyz values wrapped in a value cell */
    readonly vertexBuffer: ValueCell<Float32Array>,
    /** Group buffer as array of group ids for each vertex wrapped in a value cell */
    readonly groupBuffer: ValueCell<Float32Array>,
}

export namespace Point {
    export function createEmpty(point?: Point): Point {
        const vb = point ? point.vertexBuffer.ref.value : new Float32Array(0)
        const gb = point ? point.groupBuffer.ref.value : new Float32Array(0)
        return {
            kind: 'point',
            vertexCount: 0,
            vertexBuffer: point ? ValueCell.update(point.vertexBuffer, vb) : ValueCell.create(vb),
            groupBuffer: point ? ValueCell.update(point.groupBuffer, gb) : ValueCell.create(gb),
        }
    }

    export function transformImmediate(point: Point, t: Mat4) {
        transformRangeImmediate(point, t, 0, point.vertexCount)
    }

    export function transformRangeImmediate(point: Point, t: Mat4, offset: number, count: number) {
        const v = point.vertexBuffer.ref.value
        transformPositionArray(t, v, offset, count)
        ValueCell.update(point.vertexBuffer, v);
    }

    //

    export const DefaultProps = {
        ...Geometry.DefaultProps,
        pointSizeAttenuation: false,
        pointFilledCircle: false,
        pointEdgeBleach: 0.2,
        sizeTheme: { name: 'uniform', value: 1 } as SizeThemeProps,
    }
    export type Props = typeof DefaultProps

    export async function createValues(ctx: RuntimeContext, point: Point, transform: TransformData, locationIt: LocationIterator, props: Props): Promise<PointValues> {
        const { instanceCount, groupCount } = locationIt
        const color = await createColors(ctx, locationIt, props.colorTheme)
        const size = await createSizes(ctx, locationIt, props.sizeTheme)
        const marker = createMarkers(instanceCount * groupCount)

        const counts = { drawCount: point.vertexCount, groupCount, instanceCount }

        return {
            aPosition: point.vertexBuffer,
            aGroup: point.groupBuffer,
            ...color,
            ...size,
            ...marker,
            ...transform,

            ...Geometry.createValues(props, counts),
            dPointSizeAttenuation: ValueCell.create(props.pointSizeAttenuation),
            dPointFilledCircle: ValueCell.create(props.pointFilledCircle),
            uPointEdgeBleach: ValueCell.create(props.pointEdgeBleach),
        }
    }

    export function updateValues(values: PointValues, props: Props) {
        Geometry.updateValues(values, props)
        ValueCell.updateIfChanged(values.dPointSizeAttenuation, props.pointSizeAttenuation)
        ValueCell.updateIfChanged(values.dPointFilledCircle, props.pointFilledCircle)
        ValueCell.updateIfChanged(values.uPointEdgeBleach, props.pointEdgeBleach)
    }
}