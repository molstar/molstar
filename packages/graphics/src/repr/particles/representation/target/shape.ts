/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Geometry } from '@molstar/graphics/geo/geometry/geometry';
import type { Shape } from '@molstar/model/model/shape/shape';

/**
 * The geometry for a shape target is the shape's own (pre-built) geometry, e.g. the
 * `Mesh` produced from an OBJ file. It is instanced at each particle position.
 */
export function createShapeTargetGeometry(shape: Shape, _existing?: Geometry): Geometry {
    // Rendering accepts geometries constructed by the graphics shape factory.
    return shape.geometry as Geometry;
}
