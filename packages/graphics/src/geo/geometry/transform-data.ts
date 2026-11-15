/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { getSharedSerial } from '@molstar/core/util/array';
import { Sphere3D } from '@molstar/core/math/geometry';
import {
  calcInstanceGrid,
  createEmptyInstanceGrid,
  createPendingInstanceGrid,
  type InstanceGrid,
  moveInstanceGrid,
} from '@molstar/core/math/geometry/instance-grid';
import { InstanceTransforms, type InstancePoses } from '@molstar/core/math/geometry/instance-transforms';
import { getObjectTransforms, type InstancePosesValue } from '../../gl/transform-parts.js';
import { setTransformDataCenters } from '../../gl/renderable/util.js';

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
  /**
   * the instances as poses (a particle target's) in place of `transform` (empty then): in their
   * order, the moves composed into their `matrix` (`moveTransformData`); none with matrices. An extra transform or a
   * mirroring `matrix` turns them into matrices (`updateTransformData`)
   */
  poses: ValueCell<InstancePoses | undefined>;
  /**
   * `poses` as drawn in place of `aTransform` (empty then): the grid's order (`aInstance`), `matrix` applied
   */
  instancePoses: ValueCell<InstancePosesValue>;

  uInstanceCount: ValueCell<number>;
  instanceCount: ValueCell<number>;
  /**
   * the instance of each row: in the instances' own order a view of the shared serial array (`getSharedSerial`, read
   * only: replaced, never written), in a grid's its `cellInstance`
   */
  aInstance: ValueCell<Uint32Array>;

  hasReflection: ValueCell<boolean>;

  instanceGrid: ValueCell<InstanceGrid>;
  /**
   * the cell and batch sizes the instances' grid is built with (`calcInstanceGrid`): a cell size above 1 culls the
   * instances by cells and decimates their levels of detail by cell (`hasInstanceCells`)
   */
  instanceCellSize: ValueCell<number>;
  instanceBatchSize: ValueCell<number>;
  /** the instances' bounding sphere the box sphere at any count (`TransformValues.instanceSphereFast`) */
  instanceSphereFast?: ValueCell<boolean>;
};

/**
 * The instances' placements (`aTransform`, in draw order) by their encoding, what the CPU readers take (the bounds,
 * the position iterators, the culling of marked instances): the poses as drawn (`instancePoses`), else matrices
 */
export function getInstanceTransforms(
  t: Pick<TransformData, 'aTransform' | 'instanceCount'> & { readonly instancePoses?: ValueCell<InstancePosesValue> },
): InstanceTransforms {
  return getObjectTransforms(t);
}

const instanceRows = new WeakMap<object, { version: number; rows: Uint32Array | undefined }>();

/**
 * Each instance's row of the transforms as stored (the inverse of `aInstance`: the instance grid stores the rows in its
 * order, row `k` instance `aInstance[k]`), for what reads a transform by instance: the position iterators, whose
 * instances' values the shaders read by `aInstance` (per vertex and instance colors and sizes); undefined where the
 * rows are in the instances' order. Kept per `aInstance`'s version
 */
export function getInstanceRows(t: Pick<TransformData, 'aInstance' | 'instanceCount'>): Uint32Array | undefined {
  const cell = t.aInstance;
  const version = cell.ref.version;
  const kept = instanceRows.get(cell);
  if (kept && kept.version === version) return kept.rows;
  const instances = cell.ref.value;
  const n = t.instanceCount.ref.value;
  let identity = true;
  for (let k = 0; k < n; ++k) {
    if (instances[k] !== k) {
      identity = false;
      break;
    }
  }
  let rows: Uint32Array | undefined;
  if (!identity) {
    rows = new Uint32Array(n);
    for (let k = 0; k < n; ++k) rows[instances[k]] = k;
  }
  instanceRows.set(cell, { version, rows });
  return rows;
}

/**
 * The instances' matrices as drawn (`aTransform`'s, or built from the poses: `InstanceTransforms.getMatrix`), for the
 * readers that take matrices and run rarely (export, the debug helpers, color smoothing): a new array with poses
 */
export function getInstanceMatrices(
  t: Pick<TransformData, 'aTransform' | 'instanceCount'> & { readonly instancePoses?: ValueCell<InstancePosesValue> },
): Float32Array {
  const poses = t.instancePoses?.ref.value;
  if (poses?.kind !== 'poses') return t.aTransform.ref.value;
  return writeMatrices(poses, new Float32Array(poses.count * 16));
}

/**
 * `getInstanceMatrices` in the instances' order: matrix `i` instance `i`'s (`getInstanceRows`), as the values read by
 * instance are (`tColor`, `tSize` by `aInstance`); for the readers that pair them (the geo exports). The rows as
 * stored where they are in that order already, else a new array
 */
export function getInstanceMatricesById(
  t: Pick<TransformData, 'aTransform' | 'aInstance' | 'instanceCount'> & {
    readonly instancePoses?: ValueCell<InstancePosesValue>;
  },
): Float32Array {
  const matrices = getInstanceMatrices(t);
  const rows = getInstanceRows(t);
  if (!rows) return matrices;
  const n = t.instanceCount.ref.value;
  const out = new Float32Array(n * 16);
  for (let i = 0; i < n; ++i) out.set(matrices.subarray(rows[i] * 16, rows[i] * 16 + 16), i * 16);
  return out;
}

function writeMatrices(t: InstanceTransforms, out: Float32Array) {
  const m = Mat4();
  for (let i = 0; i < t.count; ++i) Mat4.toArray(InstanceTransforms.getMatrix(t, i, m), out, i * 16);
  return out;
}

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
    ValueCell.updateIfChanged(transformData.poses, undefined);
    ValueCell.updateIfChanged(transformData.instancePoses, undefined);

    ValueCell.update(transformData.aInstance, getSharedSerial(instanceCount));

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
      poses: ValueCell.create(undefined),
      instancePoses: ValueCell.create(undefined),
      uInstanceCount: ValueCell.create(instanceCount),
      instanceCount: ValueCell.create(instanceCount),
      aInstance: ValueCell.create(getSharedSerial(instanceCount)),
      hasReflection: ValueCell.create(hasReflection),
      instanceGrid: ValueCell.create(createEmptyInstanceGrid()),
      instanceCellSize: ValueCell.create(cellSize),
      instanceBatchSize: ValueCell.create(batchSize),
    };
  }

  updateTransformData(transformData, invariantBoundingSphere, cellSize, batchSize, oneCopy);
  return transformData;
}

/**
 * The transform data of instances placed by `poses` (in their order, without `order`; a particle target's): no
 * matrix made, the grid built from the poses (`calcInstanceGrid`'s `poses`). `transformData` updated in place where
 * given (its `matrix` kept, its extra transform cleared). `pendingGrid`: the grid not built, pending
 * (`InstanceGrid.pending`, through rebuilds and moves), the instances in their own order: for what draws them by index
 * and builds the grid where a cull needs its cells (a merged object culled on the GPU)
 */
export function createPosesTransform(
  poses: InstancePoses,
  invariantBoundingSphere: Sphere3D | undefined,
  cellSize: number,
  batchSize: number,
  transformData?: TransformData,
  pendingGrid = false,
): TransformData {
  const instanceCount = poses.count;
  if (transformData) {
    ValueCell.update(transformData.matrix, transformData.matrix.ref.value);
    ValueCell.update(transformData.transform, EmptyTransform);
    ValueCell.updateIfChanged(transformData.uInstanceCount, instanceCount);
    ValueCell.updateIfChanged(transformData.instanceCount, instanceCount);
    ValueCell.updateIfChanged(transformData.hasExtraTransform, false);
    ValueCell.update(transformData.poses, poses);
    ValueCell.update(transformData.hasReflection, false);
  } else {
    transformData = {
      aTransform: ValueCell.create(EmptyTransform),
      matrix: ValueCell.create(Mat4.identity()),
      uObjectTransform: ValueCell.create(Mat4.identity()),
      transform: ValueCell.create(EmptyTransform),
      extraTransform: ValueCell.create(new Float32Array(0)),
      hasExtraTransform: ValueCell.create(false),
      poses: ValueCell.create(poses),
      instancePoses: ValueCell.create(undefined),
      uInstanceCount: ValueCell.create(instanceCount),
      instanceCount: ValueCell.create(instanceCount),
      aInstance: ValueCell.create(getSharedSerial(instanceCount)),
      hasReflection: ValueCell.create(false),
      instanceGrid: ValueCell.create(createEmptyInstanceGrid()),
      instanceCellSize: ValueCell.create(cellSize),
      instanceBatchSize: ValueCell.create(batchSize),
    };
  }
  // pending or not as asked: `updatePosesData` keeps a pending grid pending
  const grid = transformData.instanceGrid.ref.value;
  if (pendingGrid && invariantBoundingSphere && !grid.pending) {
    ValueCell.update(
      transformData.instanceGrid,
      createPendingInstanceGrid(cellSize, batchSize, invariantBoundingSphere),
    );
  } else if (!pendingGrid && grid.pending) {
    ValueCell.update(transformData.instanceGrid, createEmptyInstanceGrid());
  }
  updateTransformData(transformData, invariantBoundingSphere, cellSize, batchSize);
  return transformData;
}

/** `transform` and `aTransform` of instances placed by poses */
const EmptyTransform = new Float32Array(0);

/**
 * The poses of `transformData` turned into matrices (`transform`, in the instances' order), where they cannot stay
 * poses: an extra transform (a matrix an instance), a mirroring `matrix` (poses do not reflect)
 */
function posesToMatrices(transformData: TransformData, poses: InstancePoses) {
  const transform = writeMatrices(poses, new Float32Array(poses.count * 16));
  ValueCell.update(transformData.poses, undefined);
  ValueCell.update(transformData.instancePoses, undefined);
  ValueCell.update(transformData.transform, transform);
  ValueCell.update(transformData.aTransform, new Float32Array(poses.count * 16));
  ValueCell.update(transformData.hasReflection, checkReflection(transform, poses.count));
}

/** `transformData`'s poses (`poses`) as drawn: the grid built from them, `matrix` composed into theirs */
function updatePosesData(
  transformData: TransformData,
  poses: InstancePoses,
  invariantBoundingSphere: Sphere3D | undefined,
  cellSize: number,
  batchSize: number,
) {
  const instanceCount = poses.count;
  const matrix = transformData.matrix.ref.value;
  const placed: InstancePoses = Mat4.isIdentity(matrix)
    ? poses
    : { ...poses, matrix: poses.matrix ? Mat4.mul(Mat4(), matrix, poses.matrix) : Mat4.copy(Mat4(), matrix) };
  const grid = transformData.instanceGrid.ref.value;
  if (hasInstanceGrid(invariantBoundingSphere, instanceCount) && grid.pending) {
    // still pending: no grid built, the instances in their own order; the grid and `aInstance` kept where they hold
    // (the poses changed alone: `Particle.withPoses`), so what is laid out by them is not laid out again
    const same =
      grid.cellSize === cellSize &&
      grid.batchSize === batchSize &&
      Sphere3D.exactEquals(grid.pending.invariantBoundingSphere, invariantBoundingSphere!) &&
      transformData.aInstance.ref.value.length === instanceCount &&
      !grid.placement;
    if (!same) {
      ValueCell.update(
        transformData.instanceGrid,
        createPendingInstanceGrid(cellSize, batchSize, invariantBoundingSphere!),
      );
      ValueCell.update(transformData.aInstance, getSharedSerial(instanceCount));
    }
    ValueCell.update(transformData.instancePoses, placed);
  } else if (hasInstanceGrid(invariantBoundingSphere, instanceCount)) {
    const aInstance = getSharedSerial(instanceCount);
    const instanceGrid = calcInstanceGrid(
      {
        instanceCount,
        instance: aInstance,
        transform: EmptyTransform,
        poses: placed,
        invariantBoundingSphere: invariantBoundingSphere!,
      },
      cellSize,
      batchSize,
    );
    ValueCell.update(transformData.instanceGrid, instanceGrid);
    ValueCell.update(transformData.aInstance, instanceGrid.cellInstance);
    ValueCell.update(transformData.instancePoses, { ...placed, order: instanceGrid.cellInstance });
  } else {
    ValueCell.update(transformData.aInstance, getSharedSerial(instanceCount));
    ValueCell.update(transformData.instancePoses, placed);
  }
  // a new value each time: what packs rows reads `aTransform`'s version
  ValueCell.update(transformData.aTransform, new Float32Array(0));
  // the grid's pass over the centers, the instances' bounding sphere's (`calculateTransformDataBoundingSphere`)
  const centers = transformData.instanceGrid.ref.value.centers;
  if (centers && hasInstanceGrid(invariantBoundingSphere, instanceCount))
    setTransformDataCenters(transformData, centers);
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
  ValueCell.update(transformData.aInstance, getSharedSerial(instanceCount));
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
  const poses = transformData.poses.ref.value;
  if (poses) {
    // no instance touched: the move composed into the poses' matrix, the grid moved
    if (Mat4.determinant(move) < 0) return false;
    const matrix = poses.matrix ? Mat4.mul(Mat4(), move, poses.matrix) : Mat4.copy(Mat4(), move);
    const moved = { ...poses, matrix };
    ValueCell.update(transformData.poses, moved);
    const grid = transformData.instanceGrid.ref.value;
    const instancePoses = transformData.instancePoses.ref.value as InstancePoses;
    // a pending grid moved too: its placement, what the GPU culler's moves read
    const moves = grid.cellCount > 0 || !!grid.pending;
    ValueCell.update(transformData.instanceGrid, moves ? moveInstanceGrid(grid, move) : grid);
    ValueCell.update(transformData.instancePoses, { ...instancePoses, matrix });
    ValueCell.update(transformData.aTransform, new Float32Array(0));
    return true;
  }
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
  ValueCell.updateIfChanged(transformData.instanceCellSize, cellSize);
  ValueCell.updateIfChanged(transformData.instanceBatchSize, batchSize);
  const poses = transformData.poses.ref.value;
  if (poses) {
    if (!transformData.hasExtraTransform.ref.value && Mat4.determinant(transformData.matrix.ref.value) >= 0) {
      updatePosesData(transformData, poses, invariantBoundingSphere, cellSize, batchSize);
      return;
    }
    posesToMatrices(transformData, poses);
  }
  const wasOneCopy = isOneCopy(transformData);
  // the rebuild reads the instances' order: the grid's one copy restored to it first (a copy for the rebuild)
  if (wasOneCopy) restoreInstanceOrder(transformData, transformData.instanceCount.ref.value);
  const keepOneCopy = oneCopy ?? wasOneCopy;
  const aTransform = transformData.aTransform.ref.value;
  const instanceCount = transformData.instanceCount.ref.value;
  // the instances in their own order (the grid's order its `cellInstance`)
  const aInstance = getSharedSerial(instanceCount);
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
    // the grid's pass over the centers, the instances' bounding sphere's (`calculateTransformDataBoundingSphere`)
    if (instanceGrid.centers) setTransformDataCenters(transformData, instanceGrid.centers);
  } else {
    ValueCell.update(transformData.aInstance, aInstance);
    ValueCell.update(transformData.aTransform, aTransform);
  }
}
