/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '../../../mol-data/util';
import { Points } from './points';

export interface PointsBuilder {
    add(x: number, y: number, z: number, group: number): void
    getPoints(): Points
}

export namespace PointsBuilder {
    export function create(initialCount = 2048, chunkSize = 1024, points?: Points): PointsBuilder {
        const centers = ChunkedArray.create(Float32Array, 3, chunkSize, points ? points.centerBuffer.ref.value : initialCount);
        const groups = ChunkedArray.create(Float32Array, 1, chunkSize, points ? points.groupBuffer.ref.value : initialCount);

        return {
            add: (x: number, y: number, z: number, group: number) => {
                ChunkedArray.add3(centers, x, y, z);
                ChunkedArray.add(groups, group);
            },
            getPoints: () => {
                const cb = ChunkedArray.compact(centers, true) as Float32Array;
                const gb = ChunkedArray.compact(groups, true) as Float32Array;
                return Points.create(cb, gb, centers.elementCount, points);
            }
        };
    }
}