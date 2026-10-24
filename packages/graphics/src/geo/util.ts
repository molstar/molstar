/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, type Mat4, type Mat3 } from '@molstar/core/math/linear-algebra';
import type { NumberArray } from '@molstar/core/util/type-helpers';
import { arrayMax } from '@molstar/core/util/array';

export function normalizeVec3Array<T extends NumberArray>(a: T, count: number) {
  for (let i = 0, il = count * 3; i < il; i += 3) {
    const x = a[i];
    const y = a[i + 1];
    const z = a[i + 2];
    const s = 1 / Math.sqrt(x * x + y * y + z * z);
    a[i] = x * s;
    a[i + 1] = y * s;
    a[i + 2] = z * s;
  }
  return a;
}

const tmpV3 = Vec3();

// The array transforms below read, transform and write inline with the matrix entries hoisted out of the loop: 1.8x
// faster than per element `Vec3.fromFloat32Array`, `Vec3.transformMat4` and `Vec3.toFloat32Array` (1M positions,
// Oct 2026), and the generic `Vec3.fromArray` / `toArray` are megamorphic app-wide. Per-element loops that call other
// code use the typed `Vec3` helpers instead.

/** In place, as `Vec3.transformMat4`; skips the w divide for an affine `t`, where w is exactly 1 */
export function transformPositionArray(t: Mat4, array: NumberArray, offset: number, count: number) {
  const m0 = t[0],
    m1 = t[1],
    m2 = t[2],
    m3 = t[3];
  const m4 = t[4],
    m5 = t[5],
    m6 = t[6],
    m7 = t[7];
  const m8 = t[8],
    m9 = t[9],
    m10 = t[10],
    m11 = t[11];
  const m12 = t[12],
    m13 = t[13],
    m14 = t[14],
    m15 = t[15];
  const end = offset + count * 3;
  // two loops: a branch on `affine` inside one is ~3x slower (V8, Oct 2026)
  if (m3 === 0 && m7 === 0 && m11 === 0 && m15 === 1) {
    for (let i = offset; i < end; i += 3) {
      const x = array[i],
        y = array[i + 1],
        z = array[i + 2];
      array[i] = m0 * x + m4 * y + m8 * z + m12;
      array[i + 1] = m1 * x + m5 * y + m9 * z + m13;
      array[i + 2] = m2 * x + m6 * y + m10 * z + m14;
    }
  } else {
    for (let i = offset; i < end; i += 3) {
      const x = array[i],
        y = array[i + 1],
        z = array[i + 2];
      const w = 1 / (m3 * x + m7 * y + m11 * z + m15 || 1.0);
      array[i] = (m0 * x + m4 * y + m8 * z + m12) * w;
      array[i + 1] = (m1 * x + m5 * y + m9 * z + m13) * w;
      array[i + 2] = (m2 * x + m6 * y + m10 * z + m14) * w;
    }
  }
}

/** In place, as `Vec3.transformMat3` */
export function transformDirectionArray(n: Mat3, array: NumberArray, offset: number, count: number) {
  const n0 = n[0],
    n1 = n[1],
    n2 = n[2];
  const n3 = n[3],
    n4 = n[4],
    n5 = n[5];
  const n6 = n[6],
    n7 = n[7],
    n8 = n[8];
  for (let i = offset, il = offset + count * 3; i < il; i += 3) {
    const x = array[i],
      y = array[i + 1],
      z = array[i + 2];
    array[i] = x * n0 + y * n3 + z * n6;
    array[i + 1] = x * n1 + y * n4 + z * n7;
    array[i + 2] = x * n2 + y * n5 + z * n8;
  }
}

/** iterate over the entire array and apply the radius to each vertex */
export function appplyRadius(vertices: NumberArray, radius: number) {
  for (let i = 0, il = vertices.length; i < il; i += 3) {
    Vec3.fromArray(tmpV3, vertices, i);
    Vec3.normalize(tmpV3, tmpV3);
    Vec3.scale(tmpV3, tmpV3, radius);
    Vec3.toArray(tmpV3, vertices, i);
  }
}

const a = Vec3();
const b = Vec3();
const c = Vec3();
const cb = Vec3();
const ab = Vec3();

/**
 * indexed vertex normals weighted by triangle areas (reads inline, as the transforms above)
 *      http://www.iquilezles.org/www/articles/normals/normals.htm
 * - normals array must contain only zeros
 */
export function computeIndexedVertexNormals<T extends NumberArray>(
  vertices: NumberArray,
  indices: NumberArray,
  normals: T,
  vertexCount: number,
  triangleCount: number,
) {
  for (let i = 0, il = triangleCount * 3; i < il; i += 3) {
    const ai = indices[i] * 3;
    const bi = indices[i + 1] * 3;
    const ci = indices[i + 2] * 3;

    const bx = vertices[bi],
      by = vertices[bi + 1],
      bz = vertices[bi + 2];
    // cb = c - b, ab = a - b, n = cb x ab, as `Vec3.sub` and `Vec3.cross`
    const cbx = vertices[ci] - bx,
      cby = vertices[ci + 1] - by,
      cbz = vertices[ci + 2] - bz;
    const abx = vertices[ai] - bx,
      aby = vertices[ai + 1] - by,
      abz = vertices[ai + 2] - bz;
    const nx = cby * abz - cbz * aby;
    const ny = cbz * abx - cbx * abz;
    const nz = cbx * aby - cby * abx;

    normals[ai] += nx;
    normals[ai + 1] += ny;
    normals[ai + 2] += nz;

    normals[bi] += nx;
    normals[bi + 1] += ny;
    normals[bi + 2] += nz;

    normals[ci] += nx;
    normals[ci + 1] += ny;
    normals[ci + 2] += nz;
  }

  return normalizeVec3Array(normals, vertexCount);
}

/**
 * vertex normals for unindexed triangle soup
 * - normals array must contain only zeros
 */
export function computeVertexNormals<T extends NumberArray>(vertices: NumberArray, normals: T, vertexCount: number) {
  for (let i = 0, il = vertexCount * 3; i < il; i += 9) {
    Vec3.fromArray(a, vertices, i);
    Vec3.fromArray(b, vertices, i + 3);
    Vec3.fromArray(c, vertices, i + 6);

    Vec3.sub(cb, c, b);
    Vec3.sub(ab, a, b);
    Vec3.cross(cb, cb, ab);

    normals[i] = cb[0];
    normals[i + 1] = cb[1];
    normals[i + 2] = cb[2];

    normals[i + 3] = cb[0];
    normals[i + 4] = cb[1];
    normals[i + 5] = cb[2];

    normals[i + 6] = cb[0];
    normals[i + 7] = cb[1];
    normals[i + 8] = cb[2];
  }

  return normalizeVec3Array(normals, vertexCount);
}

/**
 * Maps groups to data, range for group i is offsets[i] to offsets[i + 1]
 */
export type GroupMapping = {
  /** data indices */
  readonly indices: ArrayLike<number>;
  /** range for group i is offsets[i] to offsets[i + 1] */
  readonly offsets: ArrayLike<number>;
};

/**
 * The `step` parameter allows to skip over repeated values in `groups`
 */
export function createGroupMapping(groups: ArrayLike<number>, dataCount: number, step = 1): GroupMapping {
  const maxId = arrayMax(groups);
  // one bucket per group id, not per item: a shape's group ids may exceed its item count (perf-audit-2 08/B2)
  const groupCount = dataCount > 0 ? maxId + 1 : 0;

  const offsets = new Int32Array(groupCount + 1);
  const bucketFill = new Int32Array(groupCount);
  const bucketSizes = new Int32Array(groupCount);

  for (let i = 0, il = dataCount * step; i < il; i += step) ++bucketSizes[groups[i]];

  let offset = 0;
  for (let i = 0; i < groupCount; i++) {
    offsets[i] = offset;
    offset += bucketSizes[i];
  }
  offsets[groupCount] = offset;

  const indices = new Int32Array(offset);
  for (let i = 0, il = dataCount * step; i < il; i += step) {
    const g = groups[i];
    const og = offsets[g] + bucketFill[g];
    indices[og] = i;
    ++bucketFill[g];
  }

  return { indices, offsets };
}
