/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ChunkedArray } from '@molstar/core/data/util';
import { Lines } from './lines.js';
import { type Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Cage } from '../../primitive/cage.js';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const caAdd4 = ChunkedArray.add4;

const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpDir = Vec3();

export interface LinesBuilder {
  add(startX: number, startY: number, startZ: number, endX: number, endY: number, endZ: number, group: number): void;
  addVec(start: Vec3, end: Vec3, group: number): void;
  addFixedCountDashes(start: Vec3, end: Vec3, segmentCount: number, group: number): void;
  addFixedLengthDashes(start: Vec3, end: Vec3, segmentLength: number, group: number): void;
  addCage(t: Mat4, cage: Cage, group: number): void;
  getLines(): Lines;
}

/**
 * The geometry's data (`LineStride` floats per line) as 2 elements of 4 floats, in the given lines' array where
 * given; `initialCount` and `chunkSize` count lines
 */
function createLineData(initialCount: number, chunkSize: number, lines?: Lines) {
  return ChunkedArray.create(
    Float32Array,
    4,
    Math.max(1, chunkSize) * 2,
    lines ? lines.lineData.ref.value.array : Math.max(1, initialCount) * 2,
  );
}

export namespace LinesBuilder {
  /** `initialCount` and `chunkSize` count lines */
  export function create(initialCount = 2048, chunkSize = 1024, lines?: Lines): LinesBuilder {
    const data = createLineData(initialCount, chunkSize, lines);

    const add = (
      startX: number,
      startY: number,
      startZ: number,
      endX: number,
      endY: number,
      endZ: number,
      group: number,
    ) => {
      caAdd4(data, startX, startY, startZ, group);
      caAdd4(data, endX, endY, endZ, 0);
    };

    const addVec = (start: Vec3, end: Vec3, group: number) => {
      add(start[0], start[1], start[2], end[0], end[1], end[2], group);
    };

    const addFixedCountDashes = (start: Vec3, end: Vec3, segmentCount: number, group: number) => {
      const d = Vec3.distance(start, end);
      const isOdd = segmentCount % 2 !== 0;
      const s = Math.floor((segmentCount + 1) / 2);
      const step = d / (segmentCount + 0.5);

      Vec3.setMagnitude(tmpDir, Vec3.sub(tmpDir, end, start), step);
      Vec3.copy(tmpVecA, start);
      for (let j = 0; j < s; ++j) {
        Vec3.add(tmpVecA, tmpVecA, tmpDir);
        if (isOdd && j === s - 1) {
          Vec3.copy(tmpVecB, end);
        } else {
          Vec3.add(tmpVecB, tmpVecA, tmpDir);
        }
        add(tmpVecA[0], tmpVecA[1], tmpVecA[2], tmpVecB[0], tmpVecB[1], tmpVecB[2], group);
        Vec3.add(tmpVecA, tmpVecA, tmpDir);
      }
    };

    return {
      add,
      addVec,
      addFixedCountDashes,
      addFixedLengthDashes: (start: Vec3, end: Vec3, segmentLength: number, group: number) => {
        const d = Vec3.distance(start, end);
        addFixedCountDashes(start, end, d / segmentLength, group);
      },
      addCage: (t: Mat4, cage: Cage, group: number) => {
        const { vertices, edges } = cage;
        for (let i = 0, il = edges.length; i < il; i += 2) {
          Vec3.fromFloat32Array(tmpVecA, vertices, edges[i] * 3);
          Vec3.fromFloat32Array(tmpVecB, vertices, edges[i + 1] * 3);
          Vec3.transformMat4(tmpVecA, tmpVecA, t);
          Vec3.transformMat4(tmpVecB, tmpVecB, t);
          add(tmpVecA[0], tmpVecA[1], tmpVecA[2], tmpVecB[0], tmpVecB[1], tmpVecB[2], group);
        }
      },
      getLines: () => {
        const array = ChunkedArray.compactTrimmed(data) as Float32Array;
        return Lines.create(array, data.elementCount / 2, lines);
      },
    };
  }
}

//

export interface StripLinesBuilder {
  start(group: number): void;
  add(x: number, y: number, z: number): void;
  addVec(v: Vec3): void;
  end(): void;
  getLines(): Lines;
}

export namespace StripLinesBuilder {
  /**
   * Connected strips of lines: a line between consecutive points of a strip, those after the first continuing the one
   * before (`Lines.lineData`); `initialCount` and `chunkSize` count points (about the lines of the strips)
   */
  export function create(initialCount = 2048, chunkSize = 1024, lines?: Lines): StripLinesBuilder {
    const data = createLineData(initialCount, chunkSize, lines);

    let stripGroup = 0;
    let pointCount = 0;
    let prevX = 0,
      prevY = 0,
      prevZ = 0;

    const addPoint = (x: number, y: number, z: number) => {
      if (pointCount > 0) {
        caAdd4(data, prevX, prevY, prevZ, stripGroup);
        caAdd4(data, x, y, z, pointCount > 1 ? 1 : 0);
      }
      prevX = x;
      prevY = y;
      prevZ = z;
      pointCount++;
    };

    return {
      start: (group: number) => {
        stripGroup = group;
        pointCount = 0;
      },
      add: (x: number, y: number, z: number) => {
        addPoint(x, y, z);
      },
      addVec: (v: Vec3) => {
        addPoint(v[0], v[1], v[2]);
      },
      end: () => {
        pointCount = 0;
      },
      getLines: () => {
        const array = ChunkedArray.compactTrimmed(data) as Float32Array;
        return Lines.create(array, data.elementCount / 2, lines);
      },
    };
  }
}
