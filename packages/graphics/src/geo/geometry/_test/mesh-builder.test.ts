/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat3, Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { ChunkedArray } from '@molstar/core/data/util';
import { MeshBuilder } from '../mesh/mesh-builder.js';
import { addSphere } from '../mesh/builder/sphere.js';
import { addCylinder } from '../mesh/builder/cylinder.js';
import { Cylinder } from '../../primitive/cylinder.js';
import { Sphere } from '../../primitive/sphere.js';
import { BoxCage } from '../../primitive/box.js';
import { type Primitive, transformPrimitive } from '../../primitive/primitive.js';
import { transformDirectionArray, transformPositionArray } from '../../util.js';

// The inlined vertex loops must give the same bits as the Vec3 helper round trip they replace
// (`Vec3.fromArray` + `Vec3.transformMat4` / `Vec3.transformMat3` into a double `Vec3`).

function affine() {
  const m = Mat4.fromRotation(Mat4(), 0.7, Vec3.normalize(Vec3(), Vec3.create(1, -2, 0.5)));
  Mat4.scale(m, m, Vec3.create(1.3, 0.4, 2.1));
  Mat4.setTranslation(m, Vec3.create(3.1, -2, 0.5));
  return m;
}

function projective() {
  const m = affine();
  m[3] = 0.05;
  m[7] = -0.02;
  m[11] = 0.01;
  m[15] = 1.2;
  return m;
}

/** affine apart from a uniform w: must still divide */
function wScaled() {
  const m = affine();
  m[15] = 2;
  return m;
}

const matrices: [string, () => Mat4][] = [
  ['identity', () => Mat4.identity()],
  ['affine', affine],
  ['projective', projective],
  ['w-scaled', wScaled],
];

const primitives: [string, () => Primitive][] = [
  ['cylinder', () => Cylinder({ radiusTop: 0.3, radiusBottom: 0.6, radialSegments: 7, topCap: true, bottomCap: true })],
  ['sphere', () => Sphere(2)],
];

type Out = { vertices: number[]; normals: number[]; indices: number[]; groups: number[] };

function emptyOut(): Out {
  return { vertices: [], normals: [], indices: [], groups: [] };
}

/** Today's addPrimitive / addPrimitiveFlipped, through the Vec3 helpers, rounded to Float32 like the builder */
function referenceAdd(out: Out, t: Mat4, p: Primitive, group: number, flipped: boolean) {
  const offset = out.vertices.length / 3;
  const n = Mat3.directionTransform(Mat3(), t);
  const v = Vec3();
  const s = flipped ? -1 : 1;
  for (let i = 0; i < p.vertices.length; i += 3) {
    Vec3.transformMat4(v, Vec3.fromArray(v, p.vertices, i), t);
    out.vertices.push(v[0], v[1], v[2]);
    Vec3.transformMat3(v, Vec3.fromArray(v, p.normals, i), n);
    out.normals.push(s * v[0], s * v[1], s * v[2]);
    out.groups.push(group);
  }
  for (let i = 0; i < p.indices.length; i += 3) {
    const [a, b, c] = [p.indices[i], p.indices[i + 1], p.indices[i + 2]];
    if (flipped) out.indices.push(c + offset, b + offset, a + offset);
    else out.indices.push(a + offset, b + offset, c + offset);
  }
}

function stateOut(state: MeshBuilder.State): Out {
  return {
    vertices: Array.from(ChunkedArray.compact(state.vertices, true)),
    normals: Array.from(ChunkedArray.compact(state.normals, true)),
    indices: Array.from(ChunkedArray.compact(state.indices, true)),
    groups: Array.from(ChunkedArray.compact(state.groups, true)),
  };
}

/** Compares the Float32 bit patterns, so -0 and NaN count too */
function bits(a: ArrayLike<number>) {
  return Array.from(new Uint32Array(Float32Array.from(a).buffer));
}

function expectSame(actual: Out, expected: Out) {
  expect(bits(actual.vertices)).toEqual(bits(expected.vertices));
  expect(bits(actual.normals)).toEqual(bits(expected.normals));
  expect(actual.indices).toEqual(expected.indices);
  expect(actual.groups).toEqual(expected.groups);
}

describe('MeshBuilder', () => {
  for (const [mName, mat] of matrices) {
    for (const [pName, prim] of primitives) {
      for (const flipped of [false, true]) {
        it(`${flipped ? 'addPrimitiveFlipped' : 'addPrimitive'}: ${pName}, ${mName} matrix, bit for bit`, () => {
          const t = mat();
          const p = prim();
          const first = Mat4.fromTranslation(Mat4(), Vec3.create(1, 2, 3));
          const state = MeshBuilder.createState(16, 16);
          const expected = emptyOut();
          // something before, so the index offset is not zero
          state.currentGroup = 4;
          MeshBuilder.addPrimitive(state, first, p);
          referenceAdd(expected, first, p, 4, false);
          state.currentGroup = 9;
          if (flipped) MeshBuilder.addPrimitiveFlipped(state, t, p);
          else MeshBuilder.addPrimitive(state, t, p);
          referenceAdd(expected, t, p, 9, flipped);
          expectSame(stateOut(state), expected);
        });
      }
    }
  }

  for (const [mName, mat] of matrices) {
    it(`addMesh: Float32 subarrays of a mesh, ${mName} matrix, bit for bit`, () => {
      // one oversized chunk, which getMesh hands over as is
      const src = MeshBuilder.createState(1000, 1000);
      src.currentGroup = 0;
      MeshBuilder.addPrimitive(src, affine(), Sphere(1));
      src.currentGroup = 1;
      MeshBuilder.addPrimitiveFlipped(src, Mat4.identity(), Cylinder({ radialSegments: 5 }));
      const mesh = MeshBuilder.getMesh(src);
      // the buffers are longer than the counts: the subarrays must cut them
      expect(mesh.vertexBuffer.ref.value.length).toBeGreaterThan(mesh.vertexCount * 3);

      const t = mat();
      const state = MeshBuilder.createState(16, 16);
      const expected = emptyOut();
      state.currentGroup = 2;
      MeshBuilder.addPrimitive(state, Mat4.identity(), Sphere(0));
      referenceAdd(expected, Mat4.identity(), Sphere(0), 2, false);
      state.currentGroup = 3;
      MeshBuilder.addMesh(state, t, mesh);
      referenceAdd(
        expected,
        t,
        {
          vertices: mesh.vertexBuffer.ref.value.subarray(0, mesh.vertexCount * 3),
          normals: mesh.normalBuffer.ref.value.subarray(0, mesh.vertexCount * 3),
          indices: mesh.indexBuffer.ref.value.subarray(0, mesh.triangleCount * 3),
        },
        3,
        false,
      );
      expectSame(stateOut(state), expected);
    });
  }

  const ring = new Float32Array([0, 0, 0, 1, 0, 0, 1.2, 0.9, 0.1, 0.3, 1.4, -0.2, -0.8, 0.7, 0.05, -0.6, -0.5, 0]);
  const fan = [0, 1, 2, 3, 4, 5, 1];

  it('addTriangleFan and addTriangleFanWithNormal', () => {
    const state = MeshBuilder.createState(16, 16);
    state.currentGroup = 5;
    MeshBuilder.addTriangleFan(state, ring, fan);
    MeshBuilder.addTriangleFanWithNormal(state, ring, fan, Vec3.create(0, 0, 1));

    const expected = MeshBuilder.createState(16, 16);
    expected.currentGroup = 5;
    const [a, b, c] = [Vec3(), Vec3(), Vec3()];
    for (const withNormal of [false, true]) {
      Vec3.fromArray(a, ring, fan[0] * 3);
      for (let i = 2; i < fan.length; ++i) {
        Vec3.fromArray(b, ring, fan[i - 1] * 3);
        Vec3.fromArray(c, ring, fan[i] * 3);
        if (withNormal) MeshBuilder.addTriangleWithNormal(expected, a, c, b, Vec3.create(0, 0, 1));
        else MeshBuilder.addTriangle(expected, a, c, b);
      }
    }
    expectSame(stateOut(state), stateOut(expected));
  });

  it('addTriangleStrip', () => {
    const strip = [0, 1, 5, 2, 4, 3];
    const state = MeshBuilder.createState(16, 16);
    MeshBuilder.addTriangleStrip(state, ring, strip);

    const expected = MeshBuilder.createState(16, 16);
    const [a, b, c, d] = [Vec3(), Vec3(), Vec3(), Vec3()];
    Vec3.fromArray(c, ring, strip[0] * 3);
    Vec3.fromArray(d, ring, strip[1] * 3);
    for (let i = 2; i < strip.length; i += 2) {
      Vec3.copy(a, c);
      Vec3.copy(b, d);
      Vec3.fromArray(c, ring, strip[i] * 3);
      Vec3.fromArray(d, ring, strip[i + 1] * 3);
      MeshBuilder.addTriangle(expected, a, b, c);
      MeshBuilder.addTriangle(expected, b, d, c);
    }
    expectSame(stateOut(state), stateOut(expected));
  });

  for (const [mName, mat] of matrices) {
    it(`addCage: box cage, ${mName} matrix`, () => {
      const t = mat();
      const cage = BoxCage();
      const state = MeshBuilder.createState(16, 16);
      MeshBuilder.addCage(state, t, cage, 0.1, 0, 5);

      const expected = MeshBuilder.createState(16, 16);
      const [a, b] = [Vec3(), Vec3()];
      const props = { radiusTop: 0.1, radiusBottom: 0.1, radialSegments: 5 };
      for (let i = 0; i < cage.edges.length; i += 2) {
        Vec3.transformMat4(a, Vec3.fromArray(a, cage.vertices, cage.edges[i] * 3), t);
        Vec3.transformMat4(b, Vec3.fromArray(b, cage.vertices, cage.edges[i + 1] * 3), t);
        addSphere(expected, a, 0.1, 0);
        addSphere(expected, b, 0.1, 0);
        addCylinder(expected, a, b, 1, props);
      }
      expectSame(stateOut(state), stateOut(expected));
    });
  }
});

describe('transformPositionArray / transformDirectionArray / transformPrimitive', () => {
  const values = Float32Array.from({ length: 3 * 40 + 2 }, (_, i) => Math.sin(i * 1.7) * (i % 7) - 0.25 * (i % 3));

  for (const [mName, mat] of matrices) {
    it(`${mName} matrix, with an offset, bit for bit`, () => {
      const t = mat();
      const n = Mat3.directionTransform(Mat3(), t);
      const count = 39;
      const offset = 2;

      const pos = values.slice();
      transformPositionArray(t, pos, offset, count);
      const dir = values.slice();
      transformDirectionArray(n, dir, offset, count);

      const refPos = values.slice();
      const refDir = values.slice();
      const v = Vec3();
      for (let i = 0; i < count * 3; i += 3) {
        Vec3.toArray(Vec3.transformMat4(v, Vec3.fromArray(v, refPos, offset + i), t), refPos, offset + i);
        Vec3.toArray(Vec3.transformMat3(v, Vec3.fromArray(v, refDir, offset + i), n), refDir, offset + i);
      }
      expect(bits(pos)).toEqual(bits(refPos));
      expect(bits(dir)).toEqual(bits(refDir));
      // untouched outside [offset, offset + 3 * count)
      expect(bits(pos.subarray(0, offset))).toEqual(bits(values.subarray(0, offset)));
      expect(bits(pos.subarray(offset + 3 * count))).toEqual(bits(values.subarray(offset + 3 * count)));

      const p = Sphere(1);
      const refP = { vertices: p.vertices.slice(), normals: p.normals.slice() };
      transformPrimitive(p, t);
      for (let i = 0; i < refP.vertices.length; i += 3) {
        Vec3.toArray(Vec3.transformMat4(v, Vec3.fromArray(v, refP.vertices, i), t), refP.vertices, i);
        Vec3.toArray(Vec3.transformMat3(v, Vec3.fromArray(v, refP.normals, i), n), refP.normals, i);
      }
      expect(bits(p.vertices)).toEqual(bits(refP.vertices));
      expect(bits(p.normals)).toEqual(bits(refP.normals));
    });
  }
});
