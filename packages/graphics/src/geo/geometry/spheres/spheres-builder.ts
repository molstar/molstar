/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ChunkedArray } from '@molstar/core/data/util';
import { Spheres } from './spheres.js';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd4 = ChunkedArray.add4;

export interface SpheresBuilder {
  add(x: number, y: number, z: number, group: number): void;
  getSpheres(): Spheres;
}

export namespace SpheresBuilder {
  export function create(initialCount = 2048, chunkSize = 1024, spheres?: Spheres): SpheresBuilder {
    // center xyz and group per item, the geometry's array (`PositionGroupStride`), reused where given
    const positionGroup = ChunkedArray.create(
      Float32Array,
      4,
      chunkSize,
      spheres ? spheres.positionGroup.ref.value.array : initialCount,
    );

    return {
      add: (x: number, y: number, z: number, group: number) => {
        caAdd4(positionGroup, x, y, z, group);
      },
      getSpheres: () => {
        const array = ChunkedArray.compactTrimmed(positionGroup) as Float32Array;
        return Spheres.fromPositionGroup(array, positionGroup.elementCount, spheres);
      },
    };
  }
}
