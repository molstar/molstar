/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { createPrimitive, Primitive } from './primitive';
import { createCage, Cage } from './cage';

export const tetrahedronVertices: ReadonlyArray<number> = [
    0.5, 0.5, 0.5,  -0.5, -0.5, 0.5,  -0.5, 0.5, -0.5,  0.5, -0.5, -0.5
];

export const tetrahedronIndices: ReadonlyArray<number> = [
    2, 1, 0,  0, 3, 2,  1, 3, 0,  2, 3, 1
];

const tetrahedronEdges: ReadonlyArray<number> = [
    0, 1,  1, 2,  2, 0,
    0, 3,  1, 3,  2, 3,
];

let tetrahedron: Primitive;
export function Tetrahedron(): Primitive {
    if (!tetrahedron) tetrahedron = createPrimitive(tetrahedronVertices, tetrahedronIndices);
    return tetrahedron;
}

const tetrahedronCage = createCage(tetrahedronVertices, tetrahedronEdges);
export function TetrahedronCage(): Cage {
    return tetrahedronCage;
}