/**
 * Copyright (c) 2022-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { fillSerial } from '@molstar/core/util/array';
import { Box3D } from './primitives/box3d.js';
import { Vec3 } from '../linear-algebra/3d/vec3.js';
import { Mat4 } from '../linear-algebra/3d/mat4.js';
import type { PositionData } from './common.js';
import { GridLookup3D } from './lookup3d/grid.js';
import { Sphere3D } from './primitives/sphere3d.js';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4Offset = Vec3.transformMat4Offset;
const v3fromArray = Vec3.fromArray;
const b3add = Box3D.add;

type TopGrid = {
  readonly batchSize: number;
  readonly batchCount: number;
  readonly batchOffsets: Uint32Array;
  readonly batchSpheres: Float32Array;
  readonly batchCell: Uint32Array;
};

type BottomGrid = {
  readonly cellSize: number;
  readonly cellCount: number;
  readonly cellOffsets: Uint32Array;
  readonly cellSpheres: Float32Array;
  readonly cellTransform: Float32Array;
  readonly cellInstance: Float32Array;
};

export type InstanceGrid = BottomGrid &
  TopGrid & {
    /**
     * the grid this one was moved from (`moveInstanceGrid`), the first built of a line of moves: the same cells, the
     * same instance order (`cellOffsets`, `cellInstance`, the batches' cells), what was split by it still splits this
     * one; none on a built grid
     */
    readonly layout?: InstanceGrid;
    /** the move from `layout` to this grid */
    readonly placement?: Mat4;
  };

export type InstanceData = {
  instanceCount: number;
  instance: Float32Array;
  transform: Float32Array;
  invariantBoundingSphere: Sphere3D;
};

export function createEmptyInstanceGrid(): InstanceGrid {
  return {
    cellSize: 0,
    cellCount: 0,
    cellOffsets: new Uint32Array(),
    cellSpheres: new Float32Array(),
    cellTransform: new Float32Array(),
    cellInstance: new Float32Array(),

    batchSize: 0,
    batchCount: 0,
    batchOffsets: new Uint32Array(),
    batchSpheres: new Float32Array(),
    batchCell: new Uint32Array(),
  };
}

/** The offsets of one cell or batch of one instance, and the cell of one batch: read only, shared by all */
const SingleOffsets = new Uint32Array([0, 1]);
const SingleBatchCell = new Uint32Array([0]);

/**
 * The grid of a single instance, as `calcInstanceGrid` computes it (a cell and a batch, both the instance's transformed
 * sphere with the invariant radius) without the grid lookups: three arrays of its own, as objects of one instance are
 * many (a unit each)
 */
function calcSingleInstanceGrid(instanceData: InstanceData, cellSize: number, batchSize: number): InstanceGrid {
  const { instance, transform, invariantBoundingSphere } = instanceData;
  const spheres = new Float32Array(4);
  const v = v3transformMat4Offset(Vec3(), invariantBoundingSphere.center, transform, 0, 0, 0);
  spheres[0] = v[0];
  spheres[1] = v[1];
  spheres[2] = v[2];
  spheres[3] = invariantBoundingSphere.radius;
  return {
    cellSize,
    cellCount: 1,
    cellOffsets: SingleOffsets,
    cellSpheres: spheres,
    cellTransform: transform.slice(0, 16),
    cellInstance: new Float32Array([instance[0]]),

    batchSize,
    batchCount: 1,
    batchOffsets: SingleOffsets,
    batchSpheres: spheres,
    batchCell: SingleBatchCell,
  };
}

/**
 * `grid` moved by `transform` (a rigid move, or with a scale: the radii by its largest): its cells' and batches' spheres
 * and its instances' transforms (`cellTransform`) moved, its offsets and instance order shared. O(cells + instances), no
 * lookup built: what a placement moved as a whole costs instead of `calcInstanceGrid`
 */
export function moveInstanceGrid(grid: InstanceGrid, transform: Mat4): InstanceGrid {
  const scale = Mat4.getMaxScaleOnAxis(transform);
  const cellSpheres = moveSpheres(grid.cellSpheres, grid.cellCount, transform, scale);
  // a single instance's grid has one array for both
  const batchSpheres =
    grid.batchSpheres === grid.cellSpheres
      ? cellSpheres
      : moveSpheres(grid.batchSpheres, grid.batchCount, transform, scale);
  const layout = grid.layout ?? grid;
  const placement = Mat4.mul(Mat4(), transform, grid.placement ?? Mat4.identity());
  return {
    ...grid,
    cellSpheres,
    cellTransform: moveTransforms(grid.cellTransform, transform),
    batchSpheres,
    layout,
    placement,
  };
}

function moveSpheres(spheres: Float32Array, count: number, m: Mat4, scale: number) {
  const out = new Float32Array(spheres.length);
  for (let i = 0, o = 0; i < count; ++i, o += 4) {
    const x = spheres[o],
      y = spheres[o + 1],
      z = spheres[o + 2];
    out[o] = m[0] * x + m[4] * y + m[8] * z + m[12];
    out[o + 1] = m[1] * x + m[5] * y + m[9] * z + m[13];
    out[o + 2] = m[2] * x + m[6] * y + m[10] * z + m[14];
    out[o + 3] = spheres[o + 3] * scale;
  }
  return out;
}

/** `m` times each affine transform of `transforms` (column-major, 16 floats each), inlined: millions of instances */
function moveTransforms(transforms: Float32Array, m: Mat4) {
  const out = new Float32Array(transforms.length);
  const [a00, a01, a02, , a10, a11, a12, , a20, a21, a22, , a30, a31, a32] = m;
  // a column at a time (the last row of an affine transform is kept)
  for (let k = 0, kl = transforms.length; k < kl; k += 4) {
    const x = transforms[k],
      y = transforms[k + 1],
      z = transforms[k + 2],
      w = transforms[k + 3];
    out[k] = a00 * x + a10 * y + a20 * z + a30 * w;
    out[k + 1] = a01 * x + a11 * y + a21 * z + a31 * w;
    out[k + 2] = a02 * x + a12 * y + a22 * z + a32 * w;
    out[k + 3] = w;
  }
  return out;
}

export function calcInstanceGrid(instanceData: InstanceData, cellSize: number, batchSize: number): InstanceGrid {
  if (instanceData.instanceCount === 1) return calcSingleInstanceGrid(instanceData, cellSize, batchSize);
  const bottomGrid = calcBottomGrid(instanceData, cellSize);
  const topGrid = calcTopGrid(bottomGrid, batchSize);

  // reorder bottom grid so top grid has consecutive cells
  // crucial for rendering performance without base-instance support

  // console.time('calcInstanceGrid reorder');
  const cellOffsets = new Uint32Array(bottomGrid.cellOffsets.length);
  const cellSpheres = new Float32Array(bottomGrid.cellSpheres.length);
  const cellInstance = new Float32Array(bottomGrid.cellInstance.length);

  // gather transforms as float64 pairs when 8-byte aligned to halve the copy work;
  // bit-exact unless a float32 pair aliases a double NaN, which requires Inf/NaN input
  const { transform } = instanceData;
  const { cellTransform } = bottomGrid;
  const use64 = transform.byteOffset % 8 === 0 && cellTransform.byteOffset % 8 === 0;
  const src64 = use64 ? new Float64Array(transform.buffer, transform.byteOffset, transform.length >> 1) : undefined;
  const dst64 = use64
    ? new Float64Array(cellTransform.buffer, cellTransform.byteOffset, cellTransform.length >> 1)
    : undefined;

  let offset = 0;
  for (let i = 0, il = topGrid.batchCell.length; i < il; ++i) {
    const cellIdx = topGrid.batchCell[i];
    const start = bottomGrid.cellOffsets[cellIdx];
    const end = bottomGrid.cellOffsets[cellIdx + 1];
    const count = end - start;

    cellOffsets[i + 1] = cellOffsets[i] + count;
    for (let j = 0; j < 4; ++j) {
      cellSpheres[i * 4 + j] = bottomGrid.cellSpheres[cellIdx * 4 + j];
    }

    for (let j = 0; j < count; ++j) {
      const idx = start + j;
      // assumes instanceData.instance is strictly serially ordered
      const id = bottomGrid.cellInstance[idx] | 0;
      if (src64 && dst64) {
        const so = id * 8,
          dof = offset * 8;
        for (let k = 0; k < 8; ++k) {
          dst64[dof + k] = src64[so + k];
        }
      } else {
        const so = id * 16,
          dof = offset * 16;
        for (let k = 0; k < 16; ++k) {
          cellTransform[dof + k] = transform[so + k];
        }
      }
      cellInstance[offset] = id;
      offset += 1;
    }
  }
  // console.timeEnd('calcInstanceGrid reorder');

  const instanceGrid = {
    cellSize: bottomGrid.cellSize,
    cellCount: bottomGrid.cellCount,
    cellOffsets,
    cellSpheres,
    cellTransform: bottomGrid.cellTransform,
    cellInstance,

    batchSize: topGrid.batchSize,
    batchCount: topGrid.batchCount,
    batchOffsets: topGrid.batchOffsets,
    batchSpheres: topGrid.batchSpheres,
    batchCell: fillSerial(topGrid.batchCell),
  };
  // console.log(instanceGrid);

  return instanceGrid;
}

function calcBottomGrid(instanceData: InstanceData, cellSize: number): BottomGrid {
  const { instanceCount, instance, transform, invariantBoundingSphere } = instanceData;
  // console.time('calcBottomGrid grid');
  const x = new Float32Array(instanceCount);
  const y = new Float32Array(instanceCount);
  const z = new Float32Array(instanceCount);
  const indices = OrderedSet.ofBounds(0, instanceCount);

  const { center, radius } = invariantBoundingSphere;
  const rv = Vec3.create(radius, radius, radius);

  // the centers and their box inlined (`Vec3.transformMat4Offset` of affine transforms, the last row
  // 0, 0, 0, 1; `Box3D.add`): millions of instances (mesoscale particles)
  const [cx, cy, cz] = center;
  let minX = Infinity,
    minY = Infinity,
    minZ = Infinity;
  let maxX = -Infinity,
    maxY = -Infinity,
    maxZ = -Infinity;
  for (let i = 0, o = 0; i < instanceCount; ++i, o += 16) {
    const px = transform[o] * cx + transform[o + 4] * cy + transform[o + 8] * cz + transform[o + 12];
    const py = transform[o + 1] * cx + transform[o + 5] * cy + transform[o + 9] * cz + transform[o + 13];
    const pz = transform[o + 2] * cx + transform[o + 6] * cy + transform[o + 10] * cz + transform[o + 14];
    x[i] = px;
    y[i] = py;
    z[i] = pz;
    if (px < minX) minX = px;
    if (py < minY) minY = py;
    if (pz < minZ) minZ = pz;
    if (px > maxX) maxX = px;
    if (py > maxY) maxY = py;
    if (pz > maxZ) maxZ = pz;
  }
  const box = Box3D.create(Vec3.create(minX, minY, minZ), Vec3.create(maxX, maxY, maxZ));
  Box3D.expand(box, box, rv);

  const positionData: PositionData = { x, y, z, indices };
  const boundary = { box, sphere: Sphere3D.fromBox3D(Sphere3D(), box) };
  const lookup = GridLookup3D(positionData, boundary, Vec3.create(cellSize, cellSize, cellSize));
  // console.timeEnd('calcBottomGrid grid');

  const { array, offset, count } = lookup.buckets;

  const cellCount = offset.length;
  const cellOffsets = new Uint32Array(cellCount + 1);
  const cellSpheres = new Float32Array(cellCount * 4);
  // only allocated here; fully populated by the reorder in calcInstanceGrid
  const cellTransform = new Float32Array(instanceCount * 16);
  const cellInstance = new Float32Array(instanceCount);

  let k = 0;
  for (let i = 0; i < cellCount; ++i) {
    const start = offset[i];
    const size = count[i];
    cellOffsets[i] = start;
    for (let j = start, jl = start + size; j < jl; ++j) {
      cellInstance[k] = instance[array[j]];
      k += 1;
    }

    if (size === 1) {
      const idx = array[start];
      cellSpheres[i * 4] = x[idx];
      cellSpheres[i * 4 + 1] = y[idx];
      cellSpheres[i * 4 + 2] = z[idx];
      cellSpheres[i * 4 + 3] = radius;
    } else {
      // the box of the cell's centers expanded by the radius, its sphere (`Sphere3D.fromBox3D`), inlined
      let bMinX = Infinity,
        bMinY = Infinity,
        bMinZ = Infinity;
      let bMaxX = -Infinity,
        bMaxY = -Infinity,
        bMaxZ = -Infinity;
      for (let j = start, jl = start + size; j < jl; ++j) {
        const idx = array[j];
        const px = x[idx],
          py = y[idx],
          pz = z[idx];
        if (px < bMinX) bMinX = px;
        if (py < bMinY) bMinY = py;
        if (pz < bMinZ) bMinZ = pz;
        if (px > bMaxX) bMaxX = px;
        if (py > bMaxY) bMaxY = py;
        if (pz > bMaxZ) bMaxZ = pz;
      }
      const hx = (bMaxX - bMinX) / 2 + radius;
      const hy = (bMaxY - bMinY) / 2 + radius;
      const hz = (bMaxZ - bMinZ) / 2 + radius;
      cellSpheres[i * 4] = (bMinX + bMaxX) / 2;
      cellSpheres[i * 4 + 1] = (bMinY + bMaxY) / 2;
      cellSpheres[i * 4 + 2] = (bMinZ + bMaxZ) / 2;
      cellSpheres[i * 4 + 3] = Math.sqrt(hx * hx + hy * hy + hz * hz);
    }
  }
  cellOffsets[cellCount] = offset[cellCount - 1] + count[cellCount - 1];

  return {
    cellSize,
    cellCount,
    cellOffsets,
    cellSpheres,
    cellTransform,
    cellInstance,
  };
}

function calcTopGrid(bottomGrid: BottomGrid, batchSize: number): TopGrid {
  const { cellCount, cellSpheres } = bottomGrid;

  // console.time('calcTopGrid grid');
  const x = new Float32Array(cellCount);
  const y = new Float32Array(cellCount);
  const z = new Float32Array(cellCount);
  const indices = OrderedSet.ofBounds(0, cellCount);

  const box = Box3D.setEmpty(Box3D());

  const v = Vec3();
  let maxRadius = 0;
  for (let i = 0; i < cellCount; ++i) {
    const i4 = i * 4;
    v3fromArray(v, cellSpheres, i4);
    x[i] = v[0];
    y[i] = v[1];
    z[i] = v[2];
    b3add(box, v);
    maxRadius = Math.max(maxRadius, cellSpheres[i4 + 3]);
  }
  const rv = Vec3.create(maxRadius, maxRadius, maxRadius);
  Box3D.expand(box, box, rv);

  const positionData: PositionData = { x, y, z, indices };
  const boundary = { box, sphere: Sphere3D.fromBox3D(Sphere3D(), box) };
  const lookup = GridLookup3D(positionData, boundary, Vec3.create(batchSize, batchSize, batchSize));
  // console.timeEnd('calcTopGrid grid');

  const { array, offset, count } = lookup.buckets;

  const batchCount = offset.length;
  const batchOffsets = new Uint32Array(batchCount + 1);
  const batchSpheres = new Float32Array(batchCount * 4);
  const batchCell = new Uint32Array(cellCount);

  const b = Box3D();
  const s = Sphere3D();

  let k = 0;
  for (let i = 0; i < batchCount; ++i) {
    const start = offset[i];
    const size = count[i];
    batchOffsets[i] = start;
    for (let j = start, jl = start + size; j < jl; ++j) {
      batchCell[k] = array[j];
      k += 1;
    }

    if (size === 1) {
      const l = array[start];
      batchSpheres[i * 4] = cellSpheres[l * 4];
      batchSpheres[i * 4 + 1] = cellSpheres[l * 4 + 1];
      batchSpheres[i * 4 + 2] = cellSpheres[l * 4 + 2];
      batchSpheres[i * 4 + 3] = cellSpheres[l * 4 + 3];
    } else {
      Box3D.setEmpty(b);
      maxRadius = 0;
      for (let j = start, jl = start + size; j < jl; ++j) {
        const l = array[j];
        v[0] = cellSpheres[l * 4];
        v[1] = cellSpheres[l * 4 + 1];
        v[2] = cellSpheres[l * 4 + 2];
        b3add(b, v);
        maxRadius = Math.max(maxRadius, cellSpheres[l * 4 + 3]);
      }
      Vec3.set(rv, maxRadius, maxRadius, maxRadius);
      Box3D.expand(b, b, rv);
      Sphere3D.fromBox3D(s, b);
      Sphere3D.toArray(s, batchSpheres, i * 4);
    }
  }
  batchOffsets[batchCount] = offset[batchCount - 1] + count[batchCount - 1];

  return {
    batchSize,
    batchCount,
    batchOffsets,
    batchSpheres,
    batchCell,
  };
}
