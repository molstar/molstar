/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 */

import { createPrimitive, Primitive } from './primitive';
import { Cage, createCage } from './cage';

const t = (1 + Math.sqrt(5)) / 2;

const icosahedronVertices: ReadonlyArray<number> = [
    -1, t, 0, 1, t, 0, -1, -t, 0, 1, -t, 0,
    0, -1, t, 0, 1, t, 0, -1, -t, 0, 1, -t,
    t, 0, -1, t, 0, 1, -t, 0, -1, -t, 0, 1
];

const icosahedronIndices: ReadonlyArray<number> = [
    0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
    1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
    3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
    4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1,
];
/** Subset of `icosahedronIndices`, for 10 faces forming a ring around the icosahedron. */
const icosahedronRingIndices: ReadonlyArray<number> = [
    1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
    4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1,
];
/** Subset of `icosahedronIndices`, for 10 faces forming the caps of the icosahedron, complementary to `icosahedronRingIndices`. */
const icosahedronCapsIndices: ReadonlyArray<number> = [
    0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
    3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
];

// /** Subset of `icosahedronIndices` made of patches and holes of 1 or 2 faces. */
// const perforatedIcosahedronIndices1: ReadonlyArray<number> = [
//     0, 11, 5, 0, 5, 1, 0, 7, 10,
//     11, 10, 2, 7, 1, 8,
//     3, 4, 2, 3, 8, 9,
//     4, 9, 5, 6, 2, 10, 8, 6, 7,
// ];
// /** Subset of `icosahedronIndices` made of patches and holes of 1 or 2 faces, complementary to `perforatedIcosahedronIndices1`. */
// const perforatedIcosahedronIndices2: ReadonlyArray<number> = [
//     0, 1, 7, 0, 10, 11,
//     1, 5, 9, 5, 11, 4, 10, 7, 6,
//     3, 9, 4, 3, 2, 6, 3, 6, 8,
//     2, 4, 11, 9, 8, 1,
// ];

const icosahedronEdges: ReadonlyArray<number> = [
    0, 11, 5, 11, 0, 5, 1, 5, 0, 1, 1, 7, 0, 7, 7, 10, 0, 10, 10, 11,
    5, 9, 4, 11, 2, 10, 6, 7, 1, 8, 3, 9, 4, 9, 3, 4, 2, 4, 2, 3,
    2, 6, 3, 6, 6, 8, 3, 8, 8, 9, 4, 5, 2, 11, 6, 10, 7, 8, 1, 9
];

let icosahedron: Primitive;
let icosahedronRing: Primitive;
let icosahedronCaps: Primitive;
export function Icosahedron(options?: { subset: 'ring' | 'caps' | undefined }): Primitive {
    if (options?.subset === undefined) {
        return icosahedron ??= createPrimitive(icosahedronVertices, icosahedronIndices);
    } else if (options?.subset === 'ring') {
        return icosahedronRing ??= createPrimitive(icosahedronVertices, icosahedronRingIndices);
    } else {
        return icosahedronCaps ??= createPrimitive(icosahedronVertices, icosahedronCapsIndices);
    }
}

const icosahedronCage = createCage(icosahedronVertices, icosahedronEdges);
export function IcosahedronCage(): Cage {
    return icosahedronCage;
}
