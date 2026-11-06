/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Tensor, Vec3 } from '@molstar/core/math/linear-algebra';
import { Mesh } from '../mesh/mesh.js';
import { computeMarchingCubesMesh } from '../../util/marching-cubes/algorithm.js';
import { getClosureFields } from '../../util/marching-cubes/closure.js';

const N = 24;
const C = 11.5;
const R = 8;

/** a sphere (field positive inside) whose cells with x < `C` are the component's (id 0), the others ignored (-2) */
function sphereFields(interiorWithoutId: boolean) {
  const space = Tensor.Space([N, N, N], [0, 1, 2], Float32Array);
  const idSpace = Tensor.Space([N, N, N], [0, 1, 2], Int32Array);
  const field = space.create();
  const id = idSpace.create();
  for (let k = 0; k < N; ++k) {
    for (let j = 0; j < N; ++j) {
      for (let i = 0; i < N; ++i) {
        const v = R - Math.hypot(i - C, j - C, k - C);
        space.set(field, i, j, k, v);
        idSpace.set(id, i, j, k, interiorWithoutId && v > 2 ? -1 : i < C ? 0 : -2);
      }
    }
  }
  return { scalarField: Tensor.create(space, field), idField: Tensor.create(idSpace, id) };
}

async function closedHemisphere(interiorWithoutId: boolean) {
  const { scalarField, idField } = sphereFields(interiorWithoutId);
  const closure = getClosureFields(scalarField, idField, 0);
  const mesh = await computeMarchingCubesMesh({
    isoLevel: 0,
    scalarField: closure.scalarField,
    idField: closure.idField,
  }).run();
  Mesh.splitClosure(mesh, closure.isClosure!);
  return mesh;
}

function createMesh(vertices: number[], indices: number[], groups?: number[]) {
  const vertexCount = vertices.length / 3;
  const normals = new Float32Array(vertexCount * 3);
  for (let i = 0; i < vertexCount; ++i) normals[i * 3 + 2] = 1;
  const g = new Float32Array(groups ?? vertexCount);
  return Mesh.create(new Float32Array(vertices), new Uint32Array(indices), normals, g, vertexCount, indices.length / 3);
}

const BoxVertices = [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1];
const BoxIndices = [
  0, 2, 1, 0, 3, 2, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7, 4, 5, 6, 4, 6, 7,
];

function position(mesh: Mesh, v: number) {
  return Vec3.fromArray(Vec3(), mesh.vertexBuffer.ref.value, v * 3);
}

function positionKey(mesh: Mesh, v: number) {
  const p = position(mesh, v);
  return `${p[0].toFixed(4)},${p[1].toFixed(4)},${p[2].toFixed(4)}`;
}

/** directed edges by vertex position */
function directedEdgeCounts(mesh: Mesh) {
  const ib = mesh.indexBuffer.ref.value;
  const counts = new Map<string, number>();
  for (let t = 0; t < mesh.triangleCount; ++t) {
    for (let j = 0; j < 3; ++j) {
      const k = `${positionKey(mesh, ib[t * 3 + j])}|${positionKey(mesh, ib[t * 3 + ((j + 1) % 3)])}`;
      counts.set(k, (counts.get(k) ?? 0) + 1);
    }
  }
  return counts;
}

/** crossings of the ray from `origin` along `dir` with the triangles (Möller-Trumbore) */
function crossings(mesh: Mesh, origin: Vec3, dir: Vec3) {
  const ib = mesh.indexBuffer.ref.value;
  const e1 = Vec3();
  const e2 = Vec3();
  const p = Vec3();
  const s = Vec3();
  const q = Vec3();
  let count = 0;
  for (let t = 0; t < mesh.triangleCount; ++t) {
    const a = position(mesh, ib[t * 3]);
    Vec3.sub(e1, position(mesh, ib[t * 3 + 1]), a);
    Vec3.sub(e2, position(mesh, ib[t * 3 + 2]), a);
    Vec3.cross(p, dir, e2);
    const det = Vec3.dot(e1, p);
    if (Math.abs(det) < 1e-12) continue;
    Vec3.sub(s, origin, a);
    const u = Vec3.dot(s, p) / det;
    if (u < 0 || u > 1) continue;
    Vec3.cross(q, s, e1);
    const v = Vec3.dot(dir, q) / det;
    if (v < 0 || u + v > 1) continue;
    if (Vec3.dot(e2, q) / det > 0) count += 1;
  }
  return count;
}

const Directions = [
  Vec3.create(0, 0, 1),
  Vec3.create(0.013, 0.007, -1),
  Vec3.create(1, 0.013, 0.007),
  Vec3.normalize(Vec3(), Vec3.create(0.31, 0.17, 0.93)),
  Vec3.normalize(Vec3(), Vec3.create(-0.83, 0.41, -0.38)),
];

describe('closure against ignored cells', () => {
  it('leaves fields without ignored cells as they are', () => {
    const space = Tensor.Space([4, 4, 4], [0, 1, 2], Float32Array);
    const idSpace = Tensor.Space([4, 4, 4], [0, 1, 2], Int32Array);
    const scalarField = Tensor.create(space, space.create());
    const idField = Tensor.create(idSpace, idSpace.create());
    const closure = getClosureFields(scalarField, idField, 0.5);
    expect(closure.isClosure).toBeUndefined();
    expect(closure.scalarField).toBe(scalarField);
    expect(closure.idField).toBe(idField);
  });

  for (const interiorWithoutId of [false, true]) {
    it(`closes a hemisphere along the cells' border${interiorWithoutId ? ', inner cells without id' : ''}`, async () => {
      const mesh = await closedHemisphere(interiorWithoutId);
      const closureVertexCount = mesh.closureVertexCount.ref.value;
      expect(closureVertexCount).toBeGreaterThan(0);

      const counts = directedEdgeCounts(mesh);
      for (const [k, c] of counts) {
        const [a, b] = k.split('|');
        expect(c).toBe(1);
        expect(counts.get(`${b}|${a}`)).toBe(1);
      }

      const patchVertexCount = mesh.vertexCount - closureVertexCount;
      for (let v = 0; v < patchVertexCount; ++v) {
        const p = position(mesh, v);
        expect(p[0]).toBeLessThan(C + 0.5);
        expect(Math.abs(Vec3.distance(p, Vec3.create(C, C, C)) - R)).toBeLessThan(0.1);
      }
      for (let v = patchVertexCount; v < mesh.vertexCount; ++v)
        expect(Math.abs(position(mesh, v)[0] - C)).toBeLessThan(0.6);

      const ib = mesh.indexBuffer.ref.value;
      const start = Mesh.getClosureTriangleStart(ib, mesh.triangleCount, patchVertexCount);
      for (let i = 0; i < start * 3; ++i) expect(ib[i]).toBeLessThan(patchVertexCount);
      for (let i = start * 3; i < mesh.triangleCount * 3; ++i) expect(ib[i]).toBeGreaterThanOrEqual(patchVertexCount);

      for (const d of Directions) {
        expect(crossings(mesh, Vec3.create(C - 3, C + 0.3, C - 0.2), d) % 2).toBe(1);
        expect(crossings(mesh, Vec3.create(C + 3, C + 0.3, C - 0.2), d) % 2).toBe(0);
      }
    });
  }

  it('gives the closure normals pointing out of the component', async () => {
    const mesh = await closedHemisphere(false);
    const nb = mesh.normalBuffer.ref.value;
    const closureVertexCount = mesh.closureVertexCount.ref.value;
    let outward = 0;
    for (let v = mesh.vertexCount - closureVertexCount; v < mesh.vertexCount; ++v) if (nb[v * 3] > 0.5) ++outward;
    expect(outward / closureVertexCount).toBeGreaterThan(0.9);
  });
});

describe('Mesh.splitClosure', () => {
  it('resets the closure when the mesh is created again', () => {
    const mesh = createMesh(BoxVertices, BoxIndices);
    Mesh.splitClosure(mesh, (_x, _y, z) => z > 0.5);
    expect(mesh.closureVertexCount.ref.value).toBeGreaterThan(0);
    Mesh.create(
      new Float32Array(BoxVertices),
      new Uint32Array(BoxIndices),
      new Float32Array(24),
      new Float32Array(8),
      8,
      12,
      mesh,
    );
    expect(mesh.closureVertexCount.ref.value).toBe(0);
  });

  it('keeps the closure last when groups are made uniform per triangle', () => {
    const mesh = createMesh(BoxVertices, BoxIndices, [0, 1, 0, 1, 0, 1, 0, 1]);
    Mesh.splitClosure(mesh, (_x, _y, z) => z > 0.5);
    const closureVertexCount = mesh.closureVertexCount.ref.value;
    const triangleCount = mesh.triangleCount;
    const ib0 = mesh.indexBuffer.ref.value;
    const start0 = Mesh.getClosureTriangleStart(ib0, triangleCount, mesh.vertexCount - closureVertexCount);
    const closure: string[] = [];
    for (let i = start0 * 3; i < triangleCount * 3; ++i) closure.push(positionKey(mesh, ib0[i]));

    Mesh.uniformTriangleGroup(mesh);

    expect(mesh.closureVertexCount.ref.value).toBe(closureVertexCount);
    expect(mesh.triangleCount).toBeGreaterThan(triangleCount);
    const ib = mesh.indexBuffer.ref.value;
    const closureVertexStart = mesh.vertexCount - closureVertexCount;
    const start = Mesh.getClosureTriangleStart(ib, mesh.triangleCount, closureVertexStart);
    const moved: string[] = [];
    for (let i = start * 3; i < mesh.triangleCount * 3; ++i) moved.push(positionKey(mesh, ib[i]));
    expect(moved).toEqual(closure);
    for (let i = 0; i < start * 3; ++i) expect(ib[i]).toBeLessThan(closureVertexStart);

    const original = Mesh.getOriginalData(mesh)!;
    expect(original.triangleCount).toBe(triangleCount);
  });
});
