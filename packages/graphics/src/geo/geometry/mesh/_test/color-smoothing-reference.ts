/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Box3D, type Sphere3D } from '@molstar/core/math/geometry';
import { lerp } from '@molstar/core/math/interpolate';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { getVolumeTexture2dLayout } from '@molstar/graphics/repr/volume/util';

/**
 * The straightforward mesh color smoothing (as `calcMeshColorSmoothing` was before its optimization), the reference
 * the tests compare with: the grid spans the bounding-sphere box, every cell within 2 cells of a sampled vertex is
 * visited, the grid is interpolated trilinearly at every vertex. CPU only (no texture fitting), serial instance index.
 */
export interface ReferenceInput {
  vertexCount: number;
  instanceCount: number;
  groupCount: number;
  transformBuffer: Float32Array;
  positionBuffer: Float32Array;
  groupBuffer: Float32Array;
  colors: Uint8Array;
  isInstanceType: boolean;
  /** `boundingSphere` for instance types, `invariantBoundingSphere` otherwise */
  sphere: Sphere3D;
  itemSize: 1 | 4;
}

function getIndex(gridDim: Vec3, width: number, itemSize: number, x: number, y: number, z: number) {
  const [xn, yn] = gridDim;
  const column = Math.floor(((z * xn) % width) / xn);
  const row = Math.floor((z * xn) / width);
  const px = column * xn + x;
  return itemSize * (row * yn * width + y * width + px);
}

/** the reference grid, in the 2d layout of `getVolumeTexture2dLayout` */
export function referenceGrid(input: ReferenceInput, resolution: number, stride: number) {
  const { vertexCount, groupCount, positionBuffer, transformBuffer, groupBuffer, colors, isInstanceType, itemSize } =
    input;
  const box = Box3D.fromSphere3D(Box3D(), input.sphere);
  const pad = 1 + resolution;
  const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
  const scaleFactor = 1 / resolution;
  const gridDim = Box3D.size(Vec3(), Box3D.scale(Box3D(), expandedBox, scaleFactor));
  Vec3.ceil(gridDim, gridDim);
  Vec3.add(gridDim, gridDim, Vec3.create(2, 2, 2));
  const { min } = expandedBox;
  const { width, height } = getVolumeTexture2dLayout(gridDim);

  const data = new Float32Array(width * height * itemSize);
  const count = new Float32Array(width * height);
  const grid = new Uint8Array(width * height * itemSize);

  const p = 2;
  const instanceCount = isInstanceType ? input.instanceCount : 1;
  const v = Vec3();
  for (let i = 0; i < instanceCount; ++i) {
    for (let j = 0; j < vertexCount; j += stride) {
      Vec3.fromArray(v, positionBuffer, j * 3);
      if (isInstanceType) Vec3.transformMat4Offset(v, v, transformBuffer, 0, 0, i * 16);
      Vec3.scale(v, Vec3.sub(v, v, min), scaleFactor);
      const ci = (i * groupCount + groupBuffer[j]) * itemSize;
      const [x, y, z] = [Math.floor(v[0]), Math.floor(v[1]), Math.floor(v[2])];
      for (let xi = Math.max(0, x - p); xi < Math.min(gridDim[0], x + p + 2); ++xi) {
        for (let yi = Math.max(0, y - p); yi < Math.min(gridDim[1], y + p + 2); ++yi) {
          for (let zi = Math.max(0, z - p); zi < Math.min(gridDim[2], z + p + 2); ++zi) {
            const d = Math.sqrt((xi - v[0]) ** 2 + (yi - v[1]) ** 2 + (zi - v[2]) ** 2);
            if (d > p) continue;
            const s = p - d;
            const index = getIndex(gridDim, width, itemSize, xi, yi, zi);
            for (let k = 0; k < itemSize; ++k) data[index + k] += colors[ci + k] * s;
            count[index / itemSize] += s;
          }
        }
      }
    }
  }
  for (let i = 0; i < count.length; ++i) {
    for (let k = 0; k < itemSize; ++k) grid[i * itemSize + k] = Math.round(data[i * itemSize + k] / count[i]);
  }
  return { grid, gridDim, width, height, min: Vec3.clone(min), scaleFactor };
}

/** `grid` (2d layout, `width` wide) interpolated trilinearly at every vertex of every instance */
export function referenceInterpolate(
  input: ReferenceInput,
  grid: { grid: Uint8Array; gridDim: ArrayLike<number>; width: number; min: ArrayLike<number>; scaleFactor: number },
) {
  const { vertexCount, positionBuffer, transformBuffer, isInstanceType, itemSize } = input;
  const gridDim = Vec3.fromArray(Vec3(), grid.gridDim, 0);
  const instanceCount = isInstanceType ? input.instanceCount : 1;
  const out = new Uint8Array(instanceCount * vertexCount * itemSize);
  const v = Vec3();
  const v0 = Vec3();
  const v1 = Vec3();
  const vd = Vec3();
  for (let i = 0; i < instanceCount; ++i) {
    for (let j = 0; j < vertexCount; ++j) {
      Vec3.fromArray(v, positionBuffer, j * 3);
      if (isInstanceType) Vec3.transformMat4Offset(v, v, transformBuffer, 0, 0, i * 16);
      for (let k = 0; k < 3; ++k) v[k] = (v[k] - grid.min[k]) * grid.scaleFactor;
      Vec3.floor(v0, v);
      Vec3.ceil(v1, v);
      // a cell is 1 wide, or 0 on a grid plane, where the weight is 0 (not 0/0, perf-audit-2 08/B3)
      Vec3.sub(vd, v, v0);
      const at = (x: number, y: number, z: number) => getIndex(gridDim, grid.width, itemSize, x, y, z);
      for (let k = 0; k < itemSize; ++k) {
        const s = (x: number, y: number, z: number) => grid.grid[at(x, y, z) + k];
        const s00 = lerp(s(v0[0], v0[1], v0[2]), s(v1[0], v0[1], v0[2]), vd[0]);
        const s01 = lerp(s(v0[0], v0[1], v1[2]), s(v1[0], v0[1], v1[2]), vd[0]);
        const s10 = lerp(s(v0[0], v1[1], v0[2]), s(v1[0], v1[1], v0[2]), vd[0]);
        const s11 = lerp(s(v0[0], v1[1], v1[2]), s(v1[0], v1[1], v1[2]), vd[0]);
        out[(i * vertexCount + j) * itemSize + k] = lerp(lerp(s00, s10, vd[1]), lerp(s01, s11, vd[1]), vd[2]);
      }
    }
  }
  return out;
}
