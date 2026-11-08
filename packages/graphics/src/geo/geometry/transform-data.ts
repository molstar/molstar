/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { fillSerial } from '@molstar/core/util/array';
import type { Sphere3D } from '@molstar/core/math/geometry';
import {
  calcInstanceGrid,
  createEmptyInstanceGrid,
  type InstanceGrid,
  moveInstanceGrid,
} from '@molstar/core/math/geometry/instance-grid';

export type TransformData = {
  /**
   * final per-instance transform calculated for instance `i` as
   * `aTransform[i] = matrix * transform[i] * extraTransform[i]`
   */
  aTransform: ValueCell<Float32Array>;
  /** global transform, see aTransform */
  matrix: ValueCell<Mat4>;
  /**
   * transient transform of an object shown moved (a drag before its release), after `aTransform`: applied by the vertex
   * stages and the culling (the grid's spheres and the bounds transformed where tested, `getDrawnBoundingSphere`), not
   * baked, so setting it touches no instance; the identity at rest. Not a placement: what reads the instances on the CPU
   * ignores it
   */
  uObjectTransform: ValueCell<Mat4>;
  /**
   * base per-instance transform, see aTransform; kept as one copy (`createTransform`'s `oneCopy`) the grid's
   * (`cellTransform`, the same array as `aTransform`): its rows in the grid's order, row `i` instance `aInstance[i]`
   */
  transform: ValueCell<Float32Array>;
  /** additional per-instance transform, see aTransform; empty until set (64 B an instance otherwise unused) */
  extraTransform: ValueCell<Float32Array>;
  /** whether extraTransform is used, see aTransform */
  hasExtraTransform: ValueCell<boolean>;

  uInstanceCount: ValueCell<number>;
  instanceCount: ValueCell<number>;
  aInstance: ValueCell<Float32Array>;

  hasReflection: ValueCell<boolean>;

  instanceGrid: ValueCell<InstanceGrid>;
};

/** whether a transform's upper 3x3 has a negative determinant (`Mat3.determinant`, inlined: millions of instances) */
function checkReflection(t: Float32Array, instanceCount: number) {
  for (let i = 0, o = 0; i < instanceCount; i++, o += 16) {
    const a00 = t[o],
      a01 = t[o + 1],
      a02 = t[o + 2];
    const a10 = t[o + 4],
      a11 = t[o + 5],
      a12 = t[o + 6];
    const a20 = t[o + 8],
      a21 = t[o + 9],
      a22 = t[o + 10];
    const det = a00 * (a22 * a11 - a12 * a21) + a01 * (-a22 * a10 + a12 * a20) + a02 * (a21 * a10 - a11 * a20);
    if (det < 0) return true;
  }
  return false;
}

/**
 * `adopt`: `transformArray` is handed over (a new one of the caller's), kept as `transform` rather than copied; the
 * array of `transformData.transform` itself is never copied onto itself. `oneCopy`: with a grid, `transform` becomes
 * the grid's (`cellTransform`, in its order), not kept in the instances' order as well (millions of instances: a copy
 * of them less)
 */
export function createTransform(
  transformArray: Float32Array,
  instanceCount: number,
  invariantBoundingSphere: Sphere3D | undefined,
  cellSize: number,
  batchSize: number,
  transformData?: TransformData,
  adopt = false,
  oneCopy = false,
): TransformData {
  const hasReflection = checkReflection(transformArray, instanceCount);

  if (transformData) {
    ValueCell.update(transformData.matrix, transformData.matrix.ref.value);
    // not into the grid's one copy (`isOneCopy`): still read by what was built from the grid until it is rebuilt
    const current = transformData.transform.ref.value;
    const transform =
      !isOneCopy(transformData) && current.length >= instanceCount * 16
        ? current
        : adopt
          ? transformArray
          : new Float32Array(instanceCount * 16);
    if (transform !== transformArray) transform.set(transformArray);
    ValueCell.update(transformData.transform, transform);
    ValueCell.updateIfChanged(transformData.uInstanceCount, instanceCount);
    ValueCell.updateIfChanged(transformData.instanceCount, instanceCount);

    // `transform` holds the instances in their order now: not the grid's one copy any more (`isOneCopy`)
    const previous = transformData.aTransform.ref.value;
    const aTransform =
      previous !== transform && previous.length >= instanceCount * 16
        ? previous
        : new Float32Array(hasInstanceGrid(invariantBoundingSphere, instanceCount) ? 0 : instanceCount * 16);
    ValueCell.update(transformData.aTransform, aTransform);

    ValueCell.updateIfChanged(transformData.hasExtraTransform, false);

    const aInstance =
      transformData.aInstance.ref.value.length >= instanceCount
        ? transformData.aInstance.ref.value
        : new Float32Array(instanceCount);
    ValueCell.update(transformData.aInstance, fillSerial(aInstance, instanceCount));

    ValueCell.update(transformData.hasReflection, hasReflection);
  } else {
    transformData = {
      // with a grid its own array (`cellTransform`) instead: no matrix, no extra transform yet
      aTransform: ValueCell.create(
        new Float32Array(hasInstanceGrid(invariantBoundingSphere, instanceCount) ? 0 : instanceCount * 16),
      ),
      matrix: ValueCell.create(Mat4.identity()),
      uObjectTransform: ValueCell.create(Mat4.identity()),
      transform: ValueCell.create(adopt ? transformArray : new Float32Array(transformArray)),
      extraTransform: ValueCell.create(new Float32Array(0)),
      hasExtraTransform: ValueCell.create(false),
      uInstanceCount: ValueCell.create(instanceCount),
      instanceCount: ValueCell.create(instanceCount),
      aInstance: ValueCell.create(fillSerial(new Float32Array(instanceCount))),
      hasReflection: ValueCell.create(hasReflection),
      instanceGrid: ValueCell.create(createEmptyInstanceGrid()),
    };
  }

  updateTransformData(transformData, invariantBoundingSphere, cellSize, batchSize, oneCopy);
  return transformData;
}

/** whether `transformData.transform` is kept as the grid's one copy (`createTransform`'s `oneCopy`): in its order */
function isOneCopy(transformData: TransformData) {
  const transform = transformData.transform.ref.value;
  return transform.length > 0 && transform === transformData.aTransform.ref.value;
}

/**
 * `transform` in the instances' order again, from the grid's one copy (row `i` instance `aInstance[i]`), into new
 * arrays (the grid's are not changed: still read where it is): before a rebuild, which reads the instances' order
 */
function restoreInstanceOrder(transformData: TransformData, instanceCount: number) {
  const rows = transformData.transform.ref.value;
  const instances = transformData.aInstance.ref.value;
  const transform = new Float32Array(instanceCount * 16);
  for (let i = 0; i < instanceCount; ++i) {
    const s = i * 16;
    const o = instances[i] * 16;
    for (let k = 0; k < 16; ++k) transform[o + k] = rows[s + k];
  }
  ValueCell.update(transformData.transform, transform);
  ValueCell.update(transformData.aTransform, new Float32Array(instanceCount * 16));
  ValueCell.update(transformData.aInstance, fillSerial(new Float32Array(instanceCount)));
}

/**
 * The instances moved as a whole by `move` (`transform` is the instances' placement moved: a list of particles moved by
 * it), what `createTransform` of the moved transforms would give but with the grid moved (`moveInstanceGrid`), not
 * built again: its cells and instance order kept, so what was split by it moves too (the GPU culling). `transform` and
 * `aTransform` moved, `aInstance` the same. False, nothing changed, with `matrix` or `extraTransform` set: they apply
 * around `transform`, the move would not be the instances'. The bounds are the caller's
 */
export function moveTransformData(transformData: TransformData, move: Mat4): boolean {
  if (transformData.hasExtraTransform.ref.value || !Mat4.isIdentity(transformData.matrix.ref.value)) return false;
  const instanceCount = transformData.instanceCount.ref.value;
  const transform = transformData.transform.ref.value;
  const grid = transformData.instanceGrid.ref.value;
  const aTransform = transformData.aTransform.ref.value;
  if (aTransform === grid.cellTransform) {
    const moved = moveInstanceGrid(grid, move);
    ValueCell.update(transformData.instanceGrid, moved);
    ValueCell.update(transformData.aTransform, moved.cellTransform);
    // the grid's one copy (`isOneCopy`): moved with the grid, the same array again
    if (transform === aTransform) {
      ValueCell.update(transformData.transform, moved.cellTransform);
    } else {
      moveAffine(transform, instanceCount, move);
      ValueCell.update(transformData.transform, transform);
    }
  } else {
    moveAffine(transform, instanceCount, move);
    moveAffine(aTransform, instanceCount, move);
    ValueCell.update(transformData.aTransform, aTransform);
    ValueCell.update(transformData.transform, transform);
  }
  // a mirroring move flips every instance's handedness
  if (Mat4.determinant(move) < 0) {
    ValueCell.update(transformData.hasReflection, checkReflection(transformData.aTransform.ref.value, instanceCount));
  }
  return true;
}

/** `m` times each of the first `count` affine transforms of `t`, in place (`Mat4.mulOffset` inlined: millions) */
function moveAffine(t: Float32Array, count: number, m: Mat4) {
  const [a00, a01, a02, , a10, a11, a12, , a20, a21, a22, , a30, a31, a32] = m;
  for (let k = 0, kl = count * 16; k < kl; k += 4) {
    const x = t[k],
      y = t[k + 1],
      z = t[k + 2],
      w = t[k + 3];
    t[k] = a00 * x + a10 * y + a20 * z + a30 * w;
    t[k + 1] = a01 * x + a11 * y + a21 * z + a31 * w;
    t[k + 2] = a02 * x + a12 * y + a22 * z + a32 * w;
  }
}

/** whether `updateTransformData` computes an instance grid, whose `cellTransform` becomes `aTransform` */
function hasInstanceGrid(invariantBoundingSphere: Sphere3D | undefined, instanceCount: number) {
  return !!invariantBoundingSphere && instanceCount > 0;
}

const identityTransform = new Float32Array(16);
Mat4.toArray(Mat4.identity(), identityTransform, 0);

export function createIdentityTransform(transformData?: TransformData): TransformData {
  return createTransform(new Float32Array(identityTransform), 1, undefined, 0, 0, transformData);
}

/**
 * updates per-instance transform calculated for instance `i` as
 * `aTransform[i] = matrix * transform[i] * extraTransform[i]`. `oneCopy` (default: as it was, `isOneCopy`): with a
 * grid and no matrix or extra transform, `transform` the grid's `cellTransform` afterwards
 */
export function updateTransformData(
  transformData: TransformData,
  invariantBoundingSphere: Sphere3D | undefined,
  cellSize: number,
  batchSize: number,
  oneCopy?: boolean,
) {
  const wasOneCopy = isOneCopy(transformData);
  // the rebuild reads the instances' order: the grid's one copy restored to it first (a copy for the rebuild)
  if (wasOneCopy) restoreInstanceOrder(transformData, transformData.instanceCount.ref.value);
  const keepOneCopy = oneCopy ?? wasOneCopy;
  const aTransform = transformData.aTransform.ref.value;
  const aInstance = transformData.aInstance.ref.value;
  const instanceCount = transformData.instanceCount.ref.value;
  const matrix = transformData.matrix.ref.value;
  const transform = transformData.transform.ref.value;
  const extraTransform = transformData.extraTransform.ref.value;
  const hasExtraTransform = transformData.hasExtraTransform.ref.value;
  const matrixIsIdentity = Mat4.isIdentity(matrix);
  const withGrid = hasInstanceGrid(invariantBoundingSphere, instanceCount);
  // `transform` as it is: the grid reads it, its `cellTransform` the new `aTransform` (not copied first: millions of
  // instances)
  const isTransform = !hasExtraTransform && matrixIsIdentity;

  if (isTransform) {
    if (!withGrid) aTransform.set(transform.subarray(0, instanceCount * 16));
  } else if (!hasExtraTransform && !matrixIsIdentity) {
    for (let i = 0; i < instanceCount; i++) {
      const i16 = i * 16;
      Mat4.mulOffset(aTransform, matrix, transform, i16, 0, i16);
    }
  } else if (hasExtraTransform && matrixIsIdentity) {
    for (let i = 0; i < instanceCount; i++) {
      const i16 = i * 16;
      Mat4.mulOffset(aTransform, extraTransform, transform, i16, i16, i16);
    }
  } else {
    for (let i = 0; i < instanceCount; i++) {
      const i16 = i * 16;
      Mat4.mulOffset(aTransform, extraTransform, transform, i16, i16, i16);
      Mat4.mulOffset(aTransform, matrix, aTransform, i16, 0, i16);
    }
  }

  for (let i = 0; i < instanceCount; i++) {
    aInstance[i] = i;
  }

  if (withGrid) {
    const instanceGrid = calcInstanceGrid(
      {
        instanceCount,
        instance: aInstance,
        transform: isTransform ? transform : aTransform,
        invariantBoundingSphere: invariantBoundingSphere!,
      },
      cellSize,
      batchSize,
    );

    ValueCell.update(transformData.instanceGrid, instanceGrid);
    ValueCell.update(transformData.aInstance, instanceGrid.cellInstance);
    ValueCell.update(transformData.aTransform, instanceGrid.cellTransform);
    // one copy: the grid's, its rows' instances `aInstance` (the instances' order let go)
    if (keepOneCopy && isTransform) ValueCell.update(transformData.transform, instanceGrid.cellTransform);
  } else {
    ValueCell.update(transformData.aInstance, aInstance);
    ValueCell.update(transformData.aTransform, aTransform);
  }
}
