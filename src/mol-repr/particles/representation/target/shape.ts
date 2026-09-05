/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Geometry } from '../../../../mol-geo/geometry/geometry';
import { Shape } from '../../../../mol-model/shape/shape';

/**
 * The geometry for a shape target is the shape's own (pre-built) geometry, e.g. the
 * `Mesh` produced from an OBJ file. It is instanced at each particle position.
 */
export function createShapeTargetGeometry(shape: Shape, _existing?: Geometry): Geometry {
    return shape.geometry;
}
