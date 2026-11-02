/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ChunkedArray } from '@molstar/core/data/util';
import { Cylinders } from './cylinders.js';
import { Vec3 } from '@molstar/core/math/linear-algebra';

export interface CylindersBuilder {
  /**
   * @param colorMode - controls if and how theme colors are interpolated
   * - for colorMode between 0 and 1 use colorMode to interpolate
   * - for colorMode == 2 do nothing, i.e., use given theme color
   * - for colorMode == 3 use position on cylinder axis to interpolate
   */
  add(
    startX: number,
    startY: number,
    startZ: number,
    endX: number,
    endY: number,
    endZ: number,
    radiusScale: number,
    topCap: boolean,
    bottomCap: boolean,
    colorMode: number,
    group: number,
  ): void;
  addFixedCountDashes(
    start: Vec3,
    end: Vec3,
    segmentCount: number,
    radiusScale: number,
    topCap: boolean,
    bottomCap: boolean,
    stubCap: boolean,
    interpolate: boolean,
    group: number,
  ): void;
  addFixedLengthDashes(
    start: Vec3,
    end: Vec3,
    segmentLength: number,
    radiusScale: number,
    topCap: boolean,
    bottomCap: boolean,
    interpolate: boolean,
    group: number,
  ): void;
  getCylinders(): Cylinders;
}

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpDir = Vec3();

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd4 = ChunkedArray.add4;

export namespace CylindersBuilder {
  /** `initialCount` and `chunkSize` count cylinders */
  export function create(initialCount = 2048, chunkSize = 1024, cylinders?: Cylinders): CylindersBuilder {
    // the geometry's data (`CylinderStride` floats per cylinder) as 3 elements of 4 floats, reused where given
    const data = ChunkedArray.create(
      Float32Array,
      4,
      chunkSize * 3,
      cylinders ? cylinders.cylinderData.ref.value.array : initialCount * 3,
    );

    const add = (
      startX: number,
      startY: number,
      startZ: number,
      endX: number,
      endY: number,
      endZ: number,
      radiusScale: number,
      topCap: boolean,
      bottomCap: boolean,
      colorMode: number,
      group: number,
    ) => {
      const cap = (topCap ? 1 : 0) + (bottomCap ? 2 : 0);
      caAdd4(data, startX, startY, startZ, group);
      caAdd4(data, endX, endY, endZ, radiusScale);
      caAdd4(data, cap, colorMode, 0, 0);
    };

    const addFixedCountDashes = (
      start: Vec3,
      end: Vec3,
      segmentCount: number,
      radiusScale: number,
      topCap: boolean,
      bottomCap: boolean,
      stubCap: boolean,
      interpolate: boolean,
      group: number,
    ) => {
      const d = Vec3.distance(start, end);
      const isOdd = segmentCount % 2 !== 0;
      const s = Math.floor((segmentCount + 1) / 2);
      const step = d / (segmentCount + 0.5);
      let colorMode = 2.0;

      Vec3.setMagnitude(tmpDir, Vec3.sub(tmpDir, end, start), step);
      Vec3.copy(tmpVecA, start);
      for (let j = 0; j < s; ++j) {
        Vec3.add(tmpVecA, tmpVecA, tmpDir);
        if (isOdd && j === s - 1) {
          Vec3.copy(tmpVecB, end);
          if (!stubCap) bottomCap = false;
        } else {
          Vec3.add(tmpVecB, tmpVecA, tmpDir);
        }
        if (interpolate) {
          colorMode = Vec3.distance(start, tmpVecB) / (d * 2);
        }
        add(
          tmpVecA[0],
          tmpVecA[1],
          tmpVecA[2],
          tmpVecB[0],
          tmpVecB[1],
          tmpVecB[2],
          radiusScale,
          topCap,
          bottomCap,
          colorMode,
          group,
        );
        Vec3.add(tmpVecA, tmpVecA, tmpDir);
      }
    };

    return {
      add,
      addFixedCountDashes,
      addFixedLengthDashes: (
        start: Vec3,
        end: Vec3,
        segmentLength: number,
        radiusScale: number,
        topCap: boolean,
        bottomCap: boolean,
        interpolate: boolean,
        group: number,
      ) => {
        const d = Vec3.distance(start, end);
        addFixedCountDashes(start, end, d / segmentLength, radiusScale, topCap, bottomCap, true, interpolate, group);
      },
      getCylinders: () => {
        const array = ChunkedArray.compactTrimmed(data) as Float32Array;
        return Cylinders.create(array, data.elementCount / 3, cylinders);
      },
    };
  }
}
