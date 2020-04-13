/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '../primitives/box3d';
import { Sphere3D } from '../primitives/sphere3d';

export interface Result<T> {
    count: number,
    indices: T[],
    squaredDistances: number[]
}

export namespace Result {
    export function add<T>(result: Result<T>, index: T, distSq: number) {
        result.squaredDistances[result.count] = distSq;
        result.indices[result.count++] = index;
    }

    export function reset(result: Result<any>) {
        result.count = 0;
    }

    export function create<T = number>(): Result<T> {
        return { count: 0, indices: [], squaredDistances: [] };
    }

    export function copy<T = number>(out: Result<T>, result: Result<T>) {
        for (let i = 0; i < result.count; ++i) {
            out.indices[i] = result.indices[i];
            out.squaredDistances[i] = result.squaredDistances[i];
        }
        out.count = result.count;
        return out;
    }
}

export interface Lookup3D<T = number> {
    // The result is mutated with each call to find.
    find(x: number, y: number, z: number, radius: number): Result<T>,
    check(x: number, y: number, z: number, radius: number): boolean,
    readonly boundary: { readonly box: Box3D, readonly sphere: Sphere3D }
    /** transient result */
    readonly result: Result<T>
}