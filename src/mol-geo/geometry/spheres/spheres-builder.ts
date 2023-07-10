/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { Spheres } from './spheres';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd3 = ChunkedArray.add3;
const caAdd = ChunkedArray.add;

export interface SpheresBuilder {
    add(x: number, y: number, z: number, group: number): void
    getSpheres(): Spheres
}

export namespace SpheresBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, spheres?: Spheres): SpheresBuilder {
        const centers = ChunkedArray.create(Float32Array, 3, chunkSize, spheres ? spheres.centerBuffer.ref.value : initialCount);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, spheres ? spheres.groupBuffer.ref.value : initialCount);

        return {
            add: (x: number, y: number, z: number, group: number) => {
                caAdd3(centers, x, y, z);
                caAdd(groups, group);
            },
            getSpheres: () => {
                const cb = ChunkedArray.compact(centers, true) as Float32Array;
                const gb = ChunkedArray.compact(groups, true) as Float32Array;
                return Spheres.create(cb, gb, centers.elementCount, spheres);
            }
        };
    }
}