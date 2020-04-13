/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPrimitive, Primitive } from './primitive';
import { createCage, Cage } from './cage';

export const octahedronVertices: ReadonlyArray<number> = [
    0.5, 0, 0,   -0.5, 0, 0,    0, 0.5, 0,
    0, -0.5, 0,   0, 0, 0.5,    0, 0, -0.5
];

export const octahedronIndices: ReadonlyArray<number> = [
    0, 2, 4,  0, 4, 3,  0, 3, 5,
    0, 5, 2,  1, 2, 5,  1, 5, 3,
    1, 3, 4,  1, 4, 2
];

export const perforatedOctahedronIndices: ReadonlyArray<number> = [
    0, 2, 4,   0, 4, 3,
    // 0, 3, 5,   0, 5, 2,
    1, 2, 5,   1, 5, 3,
    // 1, 3, 4,   1, 4, 2
];

const octahedronEdges: ReadonlyArray<number> = [
    0, 2,  1, 3,  2, 1,  3, 0,
    0, 4,  1, 4,  2, 4,  3, 4,
    0, 5,  1, 5,  2, 5,  3, 5,
];

let octahedron: Primitive;
export function Octahedron(): Primitive {
    if (!octahedron) octahedron = createPrimitive(octahedronVertices, octahedronIndices);
    return octahedron;
}

let perforatedOctahedron: Primitive;
export function PerforatedOctahedron(): Primitive {
    if (!perforatedOctahedron) perforatedOctahedron = createPrimitive(octahedronVertices, perforatedOctahedronIndices);
    return perforatedOctahedron;
}

const octahedronCage = createCage(octahedronVertices, octahedronEdges);
export function OctahedronCage(): Cage {
    return octahedronCage;
}