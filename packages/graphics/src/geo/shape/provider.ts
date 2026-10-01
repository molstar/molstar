/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ShapeGetter } from '@molstar/graphics/repr/shape/representation';
import type { Geometry, GeometryUtils } from '../geometry/geometry.js';

export interface ShapeProvider<D, G extends Geometry, P extends Geometry.Params<G>> {
    label: string
    data: D
    params: P
    getShape: ShapeGetter<D, G, P>
    geometryUtils: GeometryUtils<G>
}
