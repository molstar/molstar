/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3, type Mat4, Mat3 } from '@molstar/core/math/linear-algebra';
import { ChunkedArray } from '@molstar/core/data/util';
import { anchorMap } from '@molstar/core/util/map-anchor';
import { Mesh } from './mesh.js';
import type { Primitive } from '../../primitive/primitive.js';
import type { Cage } from '../../primitive/cage.js';
import { addSphere } from './builder/sphere.js';
import { addCylinder } from './builder/cylinder.js';

const tmpV = Vec3();
const tmpMat3 = Mat3();
const tmpVecA = Vec3();
const tmpVecB = Vec3();
const tmpVecC = Vec3();
const tmpVecD = Vec3();

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3triangleNormal = Vec3.triangleNormal;
const v3copy = Vec3.copy;
const v3fromFloat32Array = Vec3.fromFloat32Array;
const v3transformMat4 = Vec3.transformMat4;
const mat3directionTransform = Mat3.directionTransform;
const caAdd3 = ChunkedArray.add3;
const caAdd = ChunkedArray.add;

// `Vec3.fromArray` and the `Vec3` transforms are called app-wide with many array kinds, so their keyed loads are
// megamorphic, also where V8 inlines them (2-3x slower here once a session has warmed them up). Strips, fans and cages
// read their Float32Array vertices through `Vec3.fromFloat32Array`; the primitive transform below reads and transforms
// inline, with the matrix entries hoisted out of the vertex loop.

/** Adds the primitive transformed by `t`; `flipped` negates the normals and reverses the winding order */
function addTransformedPrimitive(state: MeshBuilder.State, t: Mat4, primitive: Primitive, flipped: boolean) {
  const { vertices: va, normals: na, indices: ia } = primitive;
  const { vertices, normals, indices, groups, currentGroup } = state;
  const offset = vertices.elementCount;
  const n = mat3directionTransform(tmpMat3, t);

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
  const n0 = n[0],
    n1 = n[1],
    n2 = n[2];
  const n3 = n[3],
    n4 = n[4],
    n5 = n[5];
  const n6 = n[6],
    n7 = n[7],
    n8 = n[8];
  // w is exactly 1 for an affine matrix, so `Vec3.transformMat4`'s divide changes no bits there
  const affine = m3 === 0 && m7 === 0 && m11 === 0 && m15 === 1;

  for (let i = 0, il = va.length; i < il; i += 3) {
    // position, as `Vec3.transformMat4`
    const x = va[i],
      y = va[i + 1],
      z = va[i + 2];
    let px = m0 * x + m4 * y + m8 * z + m12;
    let py = m1 * x + m5 * y + m9 * z + m13;
    let pz = m2 * x + m6 * y + m10 * z + m14;
    if (!affine) {
      const w = 1 / (m3 * x + m7 * y + m11 * z + m15 || 1.0);
      px *= w;
      py *= w;
      pz *= w;
    }
    caAdd3(vertices, px, py, pz);
    // normal, as `Vec3.transformMat3`
    const nx = na[i],
      ny = na[i + 1],
      nz = na[i + 2];
    const qx = nx * n0 + ny * n3 + nz * n6;
    const qy = nx * n1 + ny * n4 + nz * n7;
    const qz = nx * n2 + ny * n5 + nz * n8;
    if (flipped) caAdd3(normals, -qx, -qy, -qz);
    else caAdd3(normals, qx, qy, qz);
    // group
    caAdd(groups, currentGroup);
  }
  if (flipped) {
    for (let i = 0, il = ia.length; i < il; i += 3) {
      caAdd3(indices, ia[i + 2] + offset, ia[i + 1] + offset, ia[i] + offset);
    }
  } else {
    for (let i = 0, il = ia.length; i < il; i += 3) {
      caAdd3(indices, ia[i] + offset, ia[i + 1] + offset, ia[i + 2] + offset);
    }
  }
}

export namespace MeshBuilder {
  export interface State {
    currentGroup: number;
    readonly vertices: ChunkedArray<number, 3>;
    readonly normals: ChunkedArray<number, 3>;
    readonly indices: ChunkedArray<number, 3>;
    readonly groups: ChunkedArray<number, 1>;
    readonly mesh?: Mesh;
  }

  export function createState(initialCount = 2048, chunkSize = 1024, mesh?: Mesh): State {
    return {
      currentGroup: -1,
      vertices: ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.vertexBuffer.ref.value : initialCount),
      normals: ChunkedArray.create(Float32Array, 3, chunkSize, mesh ? mesh.normalBuffer.ref.value : initialCount),
      // ~2 triangles per vertex (closed surfaces, tubes and grids); `getMesh` trims what an estimate overshoots
      indices: ChunkedArray.create(Uint32Array, 3, chunkSize * 2, mesh ? mesh.indexBuffer.ref.value : initialCount * 2),
      groups: ChunkedArray.create(Float32Array, 1, chunkSize, mesh ? mesh.groupBuffer.ref.value : initialCount),
      mesh,
    };
  }

  // keeps the state's map, which the optimized code of `addTube`, `addSheet` & co. embeds (`anchorMap`)
  anchorMap(createState(1, 1));

  export function addTriangle(state: State, a: Vec3, b: Vec3, c: Vec3) {
    const { vertices, normals, indices, groups, currentGroup } = state;
    const offset = vertices.elementCount;

    // positions
    caAdd3(vertices, a[0], a[1], a[2]);
    caAdd3(vertices, b[0], b[1], b[2]);
    caAdd3(vertices, c[0], c[1], c[2]);

    v3triangleNormal(tmpV, a, b, c);
    for (let i = 0; i < 3; ++i) {
      caAdd3(normals, tmpV[0], tmpV[1], tmpV[2]); // normal
      caAdd(groups, currentGroup); // group
    }
    caAdd3(indices, offset, offset + 1, offset + 2);
  }

  export function addTriangleWithNormal(state: State, a: Vec3, b: Vec3, c: Vec3, n: Vec3) {
    const { vertices, normals, indices, groups, currentGroup } = state;
    const offset = vertices.elementCount;

    // positions
    caAdd3(vertices, a[0], a[1], a[2]);
    caAdd3(vertices, b[0], b[1], b[2]);
    caAdd3(vertices, c[0], c[1], c[2]);

    for (let i = 0; i < 3; ++i) {
      caAdd3(normals, n[0], n[1], n[2]); // normal
      caAdd(groups, currentGroup); // group
    }
    caAdd3(indices, offset, offset + 1, offset + 2);
  }

  export function addTriangleStrip(state: State, vertices: Float32Array, indices: ArrayLike<number>) {
    v3fromFloat32Array(tmpVecC, vertices, indices[0] * 3);
    v3fromFloat32Array(tmpVecD, vertices, indices[1] * 3);
    for (let i = 2, il = indices.length; i < il; i += 2) {
      v3copy(tmpVecA, tmpVecC);
      v3copy(tmpVecB, tmpVecD);
      v3fromFloat32Array(tmpVecC, vertices, indices[i] * 3);
      v3fromFloat32Array(tmpVecD, vertices, indices[i + 1] * 3);
      addTriangle(state, tmpVecA, tmpVecB, tmpVecC);
      addTriangle(state, tmpVecB, tmpVecD, tmpVecC);
    }
  }

  export function addTriangleFan(state: State, vertices: Float32Array, indices: ArrayLike<number>) {
    v3fromFloat32Array(tmpVecA, vertices, indices[0] * 3);
    for (let i = 2, il = indices.length; i < il; ++i) {
      v3fromFloat32Array(tmpVecB, vertices, indices[i - 1] * 3);
      v3fromFloat32Array(tmpVecC, vertices, indices[i] * 3);
      addTriangle(state, tmpVecA, tmpVecC, tmpVecB);
    }
  }

  export function addTriangleFanWithNormal(
    state: State,
    vertices: Float32Array,
    indices: ArrayLike<number>,
    normal: Vec3,
  ) {
    v3fromFloat32Array(tmpVecA, vertices, indices[0] * 3);
    for (let i = 2, il = indices.length; i < il; ++i) {
      v3fromFloat32Array(tmpVecB, vertices, indices[i - 1] * 3);
      v3fromFloat32Array(tmpVecC, vertices, indices[i] * 3);
      addTriangleWithNormal(state, tmpVecA, tmpVecC, tmpVecB, normal);
    }
  }

  export function addPrimitive(state: State, t: Mat4, primitive: Primitive) {
    addTransformedPrimitive(state, t, primitive, false);
  }

  /** Flips triangle normals and winding order */
  export function addPrimitiveFlipped(state: State, t: Mat4, primitive: Primitive) {
    addTransformedPrimitive(state, t, primitive, true);
  }

  export function addCage(state: State, t: Mat4, cage: Cage, radius: number, detail: number, radialSegments: number) {
    const { vertices: va, edges: ea } = cage;
    const cylinderProps = { radiusTop: radius, radiusBottom: radius, radialSegments };
    for (let i = 0, il = ea.length; i < il; i += 2) {
      v3fromFloat32Array(tmpVecA, va, ea[i] * 3);
      v3fromFloat32Array(tmpVecB, va, ea[i + 1] * 3);
      v3transformMat4(tmpVecA, tmpVecA, t);
      v3transformMat4(tmpVecB, tmpVecB, t);
      addSphere(state, tmpVecA, radius, detail);
      addSphere(state, tmpVecB, radius, detail);
      addCylinder(state, tmpVecA, tmpVecB, 1, cylinderProps);
    }
  }

  export function addMesh(state: State, t: Mat4, mesh: Mesh) {
    addPrimitive(state, t, {
      vertices: mesh.vertexBuffer.ref.value.subarray(0, mesh.vertexCount * 3),
      normals: mesh.normalBuffer.ref.value.subarray(0, mesh.vertexCount * 3),
      indices: mesh.indexBuffer.ref.value.subarray(0, mesh.triangleCount * 3),
    });
  }

  export function getMesh(state: State): Mesh {
    const { vertices, normals, indices, groups, mesh } = state;
    const vb = ChunkedArray.compactTrimmed(vertices) as Float32Array;
    const ib = ChunkedArray.compactTrimmed(indices) as Uint32Array;
    const nb = ChunkedArray.compactTrimmed(normals) as Float32Array;
    const gb = ChunkedArray.compactTrimmed(groups) as Float32Array;
    return Mesh.create(vb, ib, nb, gb, state.vertices.elementCount, state.indices.elementCount, mesh);
  }
}
