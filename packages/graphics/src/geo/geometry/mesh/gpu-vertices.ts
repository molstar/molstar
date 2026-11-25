/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ComputeContext, StorageBuffer } from '@molstar/graphics/gl/compute';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';

/**
 * A mesh's vertices on the GPU for kernels (the smoothing's splat and sampling): its positions as `f32` xyz and its
 * groups as `u32` in storage buffers of the compute context, written once a version of its cells (`aPosition`,
 * `aGroup`), not once a layer. Not the render item's vertex buffers: the geometry pages are not storage, a slot's
 * attributes may be interleaved, and they are written at the first draw. Kept by the mesh's values until
 * `releaseGpuMeshVertices` (where the values are dropped).
 */
export interface GpuMeshVertices {
  readonly compute: ComputeContext;
  /** the positions of the cell's version `positionVersion`, `vertexCount` of them */
  readonly positions: StorageBuffer;
  /** the groups of the cell's version `groupVersion` */
  readonly groups: StorageBuffer;
  readonly vertexCount: number;
  readonly positionVersion: number;
  readonly groupVersion: number;
}

const Vertices = new WeakMap<object, GpuMeshVertices>();

/** `data` in `buffer` where it fits, else in a new buffer (`buffer` destroyed) */
function write(compute: ComputeContext, buffer: StorageBuffer | undefined, data: ArrayBufferView, label: string) {
  const byteLength = Math.max(4, data.byteLength);
  if (buffer && buffer.byteLength < byteLength) {
    buffer.destroy();
    buffer = undefined;
  }
  buffer ??= compute.createBuffer(byteLength, label);
  if (data.byteLength > 0) buffer.write(data);
  return buffer;
}

/**
 * The vertices of the mesh of `values` on the GPU of `compute`: uploaded where there are none, of another context or
 * an older version of a cell (that buffer only); a buffer kept where it fits
 */
export function getGpuMeshVertices(values: MeshValues, compute: ComputeContext): GpuMeshVertices {
  const vertexCount = values.uVertexCount.ref.value;
  const positionVersion = values.aPosition.ref.version;
  const groupVersion = values.aGroup.ref.version;
  let current = Vertices.get(values);
  // of another context (a restore, a backend switch): its buffers released with it
  if (current && current.compute !== compute) {
    releaseGpuMeshVertices(values);
    current = undefined;
  }
  const sameCount = current?.vertexCount === vertexCount;
  if (current && sameCount && current.positionVersion === positionVersion && current.groupVersion === groupVersion) {
    return current;
  }
  const positions =
    current && sameCount && current.positionVersion === positionVersion
      ? current.positions
      : write(
          compute,
          current?.positions,
          values.aPosition.ref.value.subarray(0, vertexCount * 3),
          'mesh vertices positions',
        );
  const groups =
    current && sameCount && current.groupVersion === groupVersion
      ? current.groups
      : write(compute, current?.groups, values.aGroup.ref.value.subarray(0, vertexCount), 'mesh vertices groups');
  const vertices = { compute, positions, groups, vertexCount, positionVersion, groupVersion };
  Vertices.set(values, vertices);
  return vertices;
}

/** Releases the vertices of the mesh of `values` on the GPU, if any */
export function releaseGpuMeshVertices(values: object) {
  const vertices = Vertices.get(values);
  if (!vertices) return;
  vertices.positions.destroy();
  vertices.groups.destroy();
  Vertices.delete(values);
}
