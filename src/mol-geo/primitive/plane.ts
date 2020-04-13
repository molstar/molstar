/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Primitive } from './primitive';
import { Cage } from './cage';

const plane: Primitive = {
    vertices: new Float32Array([
        -0.5, 0.5, 0,
        0.5, 0.5, 0,
        -0.5, -0.5, 0,
        0.5, -0.5, 0
    ]),
    normals: new Float32Array([
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1
    ]),
    indices: new Uint32Array([
        0, 2, 1,
        1, 2, 3
    ])
};

const planeCage: Cage = {
    vertices: plane.vertices,
    edges: new Uint32Array([ 0, 1,  2, 3,  3, 1,  2, 0 ])
};

export function Plane(): Primitive {
    return plane;
}

export function PlaneCage(): Cage {
    return planeCage;
}