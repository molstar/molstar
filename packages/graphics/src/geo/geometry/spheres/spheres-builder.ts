/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { int32ToFloatBits } from '../center-group.js';
import { ChunkedArray } from '@molstar/core/data/util';
import { Spheres } from './spheres.js';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd4 = ChunkedArray.add4;

export interface SpheresBuilder {
  add(x: number, y: number, z: number, group: number): void;
  /** `base`: the spheres every `base`-th of a whole's (`Spheres.fromPositionGroup`) */
  getSpheres(base?: number): Spheres;
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
        // the group as its bits (`rgb32f-i32`: exact past 2^24)
        caAdd4(positionGroup, x, y, z, int32ToFloatBits(group));
      },
      getSpheres: (base?: number) => {
        const array = ChunkedArray.compactTrimmed(positionGroup) as Float32Array;
        return Spheres.fromPositionGroup(array, positionGroup.elementCount, spheres, base);
      },
    };
  }
}
