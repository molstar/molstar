/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPrimitive, Primitive } from './primitive';

export const octahedronVertices: ReadonlyArray<number> = [
    0.5, 0, 0,   -0.5, 0, 0,    0, 0.5, 0,
    0, -0.5, 0,     0, 0, 0.5,  0, 0, -0.5
];

export const octahedronIndices: ReadonlyArray<number> = [
    0, 2, 4,  0, 4, 3,  0, 3, 5,
    0, 5, 2,  1, 2, 5,  1, 5, 3,
    1, 3, 4,  1, 4, 2
];

const octahedron = createPrimitive(octahedronVertices, octahedronIndices)

export function Octahedron(): Primitive { return octahedron }