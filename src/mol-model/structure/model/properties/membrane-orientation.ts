/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.bittrich@rcsb.org>
 */

import { Vec3 } from '../../../../mol-math/linear-algebra';

interface MembraneOrientation {
    // point in membrane boundary
    readonly p1: Vec3,
    // point in opposite side of membrane boundary
    readonly p2: Vec3,
    // normal vector of membrane layer
    readonly normal: Vec3,
    // the radius of the membrane layer
    readonly radius: number,
    readonly centroid: Vec3
}

function MembraneOrientation(p1: Vec3, p2: Vec3, normal: Vec3, radius: number, centroid: Vec3): MembraneOrientation {
    return { p1, p2, normal, radius, centroid };
}

export { MembraneOrientation };