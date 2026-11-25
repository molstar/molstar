/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ComputeContext, StorageBuffer } from '@molstar/graphics/gl/compute';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';

/**
 * A mesh's vertices on the GPU for kernels (the smoothing's sampling): its positions as `f32` xyz in a storage buffer of
 * the compute context, written once a version of its position cell (`aPosition`), not once a layer. Not the render
 * item's vertex buffer: the geometry pages are not storage, a slot's attributes may be interleaved, and it is written
 * at the first draw. Kept by the mesh's values until `releaseGpuMeshVertices` (where the values are dropped).
 */
export interface GpuMeshVertices {
  readonly compute: ComputeContext;
  /** the positions of the cell's version `version`, `vertexCount` of them */
  readonly positions: StorageBuffer;
  readonly vertexCount: number;
  readonly version: number;
}

const Vertices = new WeakMap<object, GpuMeshVertices>();

/**
 * The vertices of the mesh of `values` on the GPU of `compute`: uploaded where there are none, of another context or
 * an older version of the positions; the buffer kept where it fits
 */
export function getGpuMeshVertices(values: MeshValues, compute: ComputeContext): GpuMeshVertices {
  const cell = values.aPosition;
  const vertexCount = values.uVertexCount.ref.value;
  const current = Vertices.get(values);
  if (current?.compute === compute && current.version === cell.ref.version && current.vertexCount === vertexCount) {
    return current;
  }
  const data = cell.ref.value.subarray(0, vertexCount * 3);
  const byteLength = Math.max(4, data.byteLength);
  let positions = current?.compute === compute ? current.positions : undefined;
  if (positions && positions.byteLength < byteLength) {
    positions.destroy();
    positions = undefined;
  }
  // of another context (a restore, a backend switch): its buffer released with it
  if (current && current.compute !== compute) current.positions.destroy();
  positions ??= compute.createBuffer(byteLength, 'mesh vertices positions');
  if (data.byteLength > 0) positions.write(data);
  const vertices = { compute, positions, vertexCount, version: cell.ref.version };
  Vertices.set(values, vertices);
  return vertices;
}

/** Releases the vertices of the mesh of `values` on the GPU, if any */
export function releaseGpuMeshVertices(values: object) {
  const vertices = Vertices.get(values);
  if (!vertices) return;
  vertices.positions.destroy();
  Vertices.delete(values);
}
