/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPrimitive, Primitive } from './primitive';
import { createCage, Cage } from './cage';

export const tetrahedronVertices: ReadonlyArray<number> = [
    0.7071, 0, 0,  -0.3535, 0.6123, 0,  -0.3535, -0.6123, 0,
    0, 0, 0.7071,  0, 0, -0.7071

];

export const tetrahedronIndices: ReadonlyArray<number> = [
    4, 1, 0,  4, 2, 1,  4, 0, 2,
    0, 1, 3,  1, 2, 3,  2, 0, 3,
];

const tetrahedronEdges: ReadonlyArray<number> = [
    0, 1,  1, 2,  2, 0,
    0, 3,  1, 3,  2, 3,
    0, 4,  1, 4,  2, 4,
]

let tetrahedron: Primitive
export function Tetrahedron(): Primitive {
    if (!tetrahedron) tetrahedron = createPrimitive(tetrahedronVertices, tetrahedronIndices)
    return tetrahedron
}

const tetrahedronCage = createCage(tetrahedronVertices, tetrahedronEdges)
export function TetrahedronCage(): Cage {
    return tetrahedronCage
}