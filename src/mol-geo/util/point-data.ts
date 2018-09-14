/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from 'mol-util/value-cell'
import { Point } from '../geometry/point/point';

type PointData = {
    aPosition: ValueCell<Float32Array>,
    aGroup: ValueCell<Float32Array>,
}

export function getPointData(point: Point): PointData {
    return {
        aPosition: point.vertexBuffer,
        aGroup: point.groupBuffer,
    }
}