/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mesh } from '../mesh/mesh.js';

function createMesh(vertices: number[], indices: number[]) {
  const vertexCount = vertices.length / 3;
  const normals = new Float32Array(vertexCount * 3);
  for (let i = 0; i < vertexCount; ++i) normals[i * 3 + 2] = 1;
  return Mesh.create(
    new Float32Array(vertices),
    new Uint32Array(indices),
    normals,
    new Float32Array(vertexCount),
    vertexCount,
    indices.length / 3,
  );
}

describe('Mesh.smoothEdges', () => {
  it('keeps border vertices finite where border neighbours coincide', () => {
    const vertices = [0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0.5, 0.5, 0];
    const mesh = createMesh(vertices, [5, 0, 1, 5, 1, 2, 5, 2, 3, 5, 3, 4, 5, 4, 0]);
    Mesh.smoothEdges(mesh, { iterations: 2, maxNewEdgeLength: Math.sqrt(2) });
    const vb = mesh.vertexBuffer.ref.value;
    for (let i = 0; i < mesh.vertexCount * 3; ++i) expect(Number.isFinite(vb[i])).toBe(true);
  });
});
