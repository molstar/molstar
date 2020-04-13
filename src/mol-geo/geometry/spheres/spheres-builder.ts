/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { Spheres } from './spheres';

const quadMapping = new Float32Array([
    -1.0, 1.0,
    -1.0, -1.0,
    1.0, 1.0,
    1.0, -1.0
]);

const quadIndices = new Uint16Array([
    0, 1, 2,
    1, 3, 2
]);

export interface SpheresBuilder {
    add(x: number, y: number, z: number, group: number): void
    getSpheres(): Spheres
}

export namespace SpheresBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, spheres?: Spheres): SpheresBuilder {
        initialCount *= 4;
        chunkSize *= 4;
        const centers = ChunkedArray.create(Float32Array, 3, chunkSize, spheres ? spheres.centerBuffer.ref.value : initialCount);
        const mappings = ChunkedArray.create(Float32Array, 2, chunkSize, spheres ? spheres.mappingBuffer.ref.value : initialCount);
        const indices = ChunkedArray.create(Uint32Array, 3, chunkSize / 2, spheres ? spheres.indexBuffer.ref.value : initialCount / 2);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, spheres ? spheres.groupBuffer.ref.value : initialCount);

        return {
            add: (x: number, y: number, z: number, group: number) => {
                const offset = centers.elementCount;
                for (let i = 0; i < 4; ++i) {
                    ChunkedArray.add3(centers, x, y, z);
                    ChunkedArray.add2(mappings, quadMapping[i * 2], quadMapping[i * 2 + 1]);
                    ChunkedArray.add(groups, group);
                }
                ChunkedArray.add3(indices, offset + quadIndices[0], offset + quadIndices[1], offset + quadIndices[2]);
                ChunkedArray.add3(indices, offset + quadIndices[3], offset + quadIndices[4], offset + quadIndices[5]);
            },
            getSpheres: () => {
                const cb = ChunkedArray.compact(centers, true) as Float32Array;
                const mb = ChunkedArray.compact(mappings, true) as Float32Array;
                const ib = ChunkedArray.compact(indices, true) as Uint32Array;
                const gb = ChunkedArray.compact(groups, true) as Float32Array;
                return Spheres.create(cb, mb, ib, gb, centers.elementCount / 4, spheres);
            }
        };
    }
}