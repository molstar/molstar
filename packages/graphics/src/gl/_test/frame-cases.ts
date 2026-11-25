/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The objects the allocation cases draw (frame.quick.ts): small meshes and sphere sets, each a draw.

import { ColorNames } from '@molstar/core/util/color/names';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mesh } from '../../geo/geometry/mesh/mesh.js';
import { Spheres } from '../../geo/geometry/spheres/spheres.js';
import { createRenderObject } from '../render-object.js';

/** A mesh of `vertexCount` vertices in triangles along x from `offset` */
export function createMeshObject(vertexCount: number, offset: number) {
  const vertices = new Float32Array(vertexCount * 3);
  const normals = new Float32Array(vertexCount * 3);
  const groups = new Float32Array(vertexCount);
  for (let i = 0; i < vertexCount; ++i) {
    vertices[i * 3] = (i % 3) + offset;
    vertices[i * 3 + 1] = Math.floor(i / 3) % 2;
    normals[i * 3 + 2] = 1;
  }
  const triangleCount = Math.floor(vertexCount / 3);
  const indices = new Uint32Array(triangleCount * 3);
  for (let i = 0, il = triangleCount * 3; i < il; ++i) indices[i] = i;
  const mesh = Mesh.create(vertices, indices, normals, groups, vertexCount, triangleCount);
  const props = PD.getDefaultValues(Mesh.Params);
  const values = Mesh.Utils.createValuesSimple(mesh, props, ColorNames.orange, 1);
  const state = Mesh.Utils.createRenderableState(props);
  return createRenderObject('mesh', values, state, -1);
}

/** `count` spheres along y from `offset` */
export function createSpheresObject(count: number, offset: number) {
  const centers = new Float32Array(count * 3);
  const groups = new Float32Array(count);
  for (let i = 0; i < count; ++i) {
    centers[i * 3] = offset;
    centers[i * 3 + 1] = i;
    groups[i] = i;
  }
  const spheres = Spheres.create(centers, groups, count);
  const props = PD.getDefaultValues(Spheres.Params);
  const values = Spheres.Utils.createValuesSimple(spheres, props, ColorNames.orange, 1);
  const state = Spheres.Utils.createRenderableState(props);
  return createRenderObject('spheres', values, state, -1);
}
