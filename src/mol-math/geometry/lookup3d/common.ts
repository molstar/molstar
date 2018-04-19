/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Box3D } from '../primitives/box3d'
import { Sphere3D } from '../primitives/sphere3d'

export interface Result {
    count: number,
    indices: number[],
    squaredDistances: number[]
}

export namespace Result {
    export function add(result: Result, index: number, distSq: number) {
        result.squaredDistances[result.count] = distSq;
        result.indices[result.count++] = index;
    }

    export function reset(result: Result) {
        result.count = 0;
    }

    export function create(): Result {
        return { count: 0, indices: [], squaredDistances: [] };
    }
}

export interface Lookup3D {
    // The result is mutated with each call to find.
    find(x: number, y: number, z: number, radius: number): Result,
    check(x: number, y: number, z: number, radius: number): boolean,
    boundingBox: Box3D,
    boundingSphere: Sphere3D
}