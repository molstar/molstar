/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { fillIdentityTransform } from '@molstar/core/math/linear-algebra/3d/transform-array';
import { ValueCell } from '@molstar/core/util';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Sphere3D } from '@molstar/core/math/geometry';
import {
  createIdentityTransform,
  createPosesTransform,
  createTransform,
  getInstanceMatrices,
  moveTransformData,
  type TransformData,
  updateTransformData,
} from '../transform-data.js';
import { InstanceTransforms, type InstancePoses } from '@molstar/core/math/geometry/instance-transforms';
import { Quat } from '@molstar/core/math/linear-algebra/3d/quat';
import { hasInstanceCells, writeAffineRowsOf } from '../../../gl/transform-parts.js';

function makeTransforms(count: number, factory: (i: number) => Mat4) {
  const array = new Float32Array(count * 16);
  for (let i = 0; i < count; i++) {
    Mat4.toArray(factory(i), array, i * 16);
  }
  return array;
}

function getInstance(array: Float32Array, i: number) {
  return Mat4.fromArray(Mat4(), array, i * 16);
}

describe('transform-data', () => {
  describe('fillIdentityTransform', () => {
    it('fills every instance slot with the identity matrix', () => {
      const out = fillIdentityTransform(new Float32Array(3 * 16), 3);
      for (let i = 0; i < 3; i++) {
        expect(Mat4.areEqual(getInstance(out, i), Mat4.identity(), 1e-6)).toBe(true);
      }
    });
  });

  describe('createTransform', () => {
    it('creates a TransformData with identity matrix and no extra transform', () => {
      const transformArray = makeTransforms(2, (i) => Mat4.fromTranslation(Mat4(), Vec3.create(i + 1, 0, 0)));
      const td = createTransform(transformArray, 2, undefined, 0, 0);

      expect(td.instanceCount.ref.value).toBe(2);
      expect(td.uInstanceCount.ref.value).toBe(2);
      expect(td.hasExtraTransform.ref.value).toBe(false);
      expect(Mat4.isIdentity(td.matrix.ref.value)).toBe(true);
      expect(td.hasReflection.ref.value).toBe(false);

      // matrix is identity and no extra transform => aTransform is a direct copy of transform
      expect(Array.from(td.aTransform.ref.value)).toEqual(Array.from(transformArray));
    });

    it('detects reflection (negative determinant) instance transforms', () => {
      const reflected = Mat4.fromScaling(Mat4(), Vec3.create(-1, 1, 1));
      const transformArray = makeTransforms(1, () => reflected);
      const td = createTransform(transformArray, 1, undefined, 0, 0);
      expect(td.hasReflection.ref.value).toBe(true);
    });

    it('allocates no extra transforms until they are set', () => {
      const td = createTransform(
        makeTransforms(3, () => Mat4.identity()),
        3,
        undefined,
        0,
        0,
      );
      expect(td.extraTransform.ref.value.length).toBe(0);
      expect(td.hasExtraTransform.ref.value).toBe(false);
    });

    it('resets hasExtraTransform to false when reusing an existing TransformData', () => {
      const transformArray = makeTransforms(1, () => Mat4.identity());
      const td = createTransform(transformArray, 1, undefined, 0, 0);

      ValueCell.update(td.extraTransform, fillIdentityTransform(new Float32Array(16), 1));
      ValueCell.updateIfChanged(td.hasExtraTransform, true);
      expect(td.hasExtraTransform.ref.value).toBe(true);

      createTransform(transformArray, 1, undefined, 0, 0, td);
      expect(td.hasExtraTransform.ref.value).toBe(false);
    });

    it('adopts the given array on an update that grows past the current one, copies it otherwise', () => {
      const translated = (i: number) => Mat4.fromTranslation(Mat4(), Vec3.create(i, 0, 0));
      const td = createTransform(makeTransforms(2, translated), 2, undefined, 0, 0);
      // more instances than it holds: the given array kept, not copied into a new one
      const grown = makeTransforms(4, translated);
      createTransform(grown, 4, undefined, 0, 0, td, true);
      expect(td.transform.ref.value).toBe(grown);
      // not adopted: a copy of its own
      const copied = makeTransforms(8, translated);
      createTransform(copied, 8, undefined, 0, 0, td);
      expect(td.transform.ref.value).not.toBe(copied);
      expect(getInstance(td.transform.ref.value, 7)).toEqual(translated(7));
    });
  });

  describe('createIdentityTransform', () => {
    it('creates a single identity instance transform', () => {
      const td = createIdentityTransform();
      expect(td.instanceCount.ref.value).toBe(1);
      expect(td.hasReflection.ref.value).toBe(false);
      expect(Mat4.areEqual(getInstance(td.aTransform.ref.value, 0), Mat4.identity(), 1e-6)).toBe(true);
    });
  });

  describe('updateTransformData', () => {
    const instanceCount = 3;

    function setup() {
      const transformArray = makeTransforms(instanceCount, (i) =>
        Mat4.fromTranslation(Mat4(), Vec3.create(i + 1, 0, 0)),
      );
      const td = createTransform(transformArray, instanceCount, undefined, 0, 0);
      return td;
    }

    it('copies transform directly when matrix is identity and there is no extra transform', () => {
      const td = setup();
      updateTransformData(td, undefined, 0, 0);

      for (let i = 0; i < instanceCount; i++) {
        const expected = getInstance(td.transform.ref.value, i);
        expect(Mat4.areEqual(getInstance(td.aTransform.ref.value, i), expected, 1e-6)).toBe(true);
      }
    });

    it('applies the global matrix per-instance when there is no extra transform', () => {
      const td = setup();
      const matrix = Mat4.fromScaling(Mat4(), Vec3.create(2, 3, 4));
      Mat4.copy(td.matrix.ref.value, matrix);
      ValueCell.update(td.matrix, td.matrix.ref.value);

      updateTransformData(td, undefined, 0, 0);

      for (let i = 0; i < instanceCount; i++) {
        const transform_i = getInstance(td.transform.ref.value, i);
        const expected = Mat4.mul(Mat4(), matrix, transform_i);
        const actual = getInstance(td.aTransform.ref.value, i);
        expect(Mat4.areEqual(actual, expected, 1e-6)).toBe(true);
      }
    });

    it('applies per-instance extra transform when matrix is identity', () => {
      const td = setup();
      const extraTransformArray = makeTransforms(instanceCount, (i) =>
        Mat4.fromTranslation(Mat4(), Vec3.create(0, i + 1, 0)),
      );
      ValueCell.update(td.extraTransform, extraTransformArray);
      ValueCell.updateIfChanged(td.hasExtraTransform, true);

      updateTransformData(td, undefined, 0, 0);

      for (let i = 0; i < instanceCount; i++) {
        const transform_i = getInstance(td.transform.ref.value, i);
        const extra_i = getInstance(extraTransformArray, i);
        const expected = Mat4.mul(Mat4(), extra_i, transform_i);
        const actual = getInstance(td.aTransform.ref.value, i);
        expect(Mat4.areEqual(actual, expected, 1e-6)).toBe(true);
      }
    });

    it('combines global matrix and per-instance extra transform', () => {
      const td = setup();
      const extraTransformArray = makeTransforms(instanceCount, (i) =>
        Mat4.fromTranslation(Mat4(), Vec3.create(0, i + 1, 0)),
      );
      ValueCell.update(td.extraTransform, extraTransformArray);
      ValueCell.updateIfChanged(td.hasExtraTransform, true);

      const matrix = Mat4.fromScaling(Mat4(), Vec3.create(2, 3, 4));
      Mat4.copy(td.matrix.ref.value, matrix);
      ValueCell.update(td.matrix, td.matrix.ref.value);

      updateTransformData(td, undefined, 0, 0);

      for (let i = 0; i < instanceCount; i++) {
        const transform_i = getInstance(td.transform.ref.value, i);
        const extra_i = getInstance(extraTransformArray, i);
        const expected = Mat4.mul(Mat4(), matrix, Mat4.mul(Mat4(), extra_i, transform_i));
        const actual = getInstance(td.aTransform.ref.value, i);
        expect(Mat4.areEqual(actual, expected, 1e-6)).toBe(true);
      }
    });

    it('fills aInstance with serial indices', () => {
      const td = setup();
      updateTransformData(td, undefined, 0, 0);
      expect(Array.from(td.aInstance.ref.value.slice(0, instanceCount))).toEqual([0, 1, 2]);
    });
  });
});

function transforms(count: number) {
  const out = new Float32Array(count * 16);
  for (let i = 0; i < count; ++i) {
    Mat4.toArray(
      Mat4.fromTranslation(Mat4(), Vec3.create((i % 7) * 9, Math.floor(i / 7) * 5, (i % 3) * 4)),
      out,
      i * 16,
    );
  }
  return out;
}

const move = Mat4.mul(
  Mat4(),
  Mat4.fromTranslation(Mat4(), Vec3.create(20, -3, 7)),
  Mat4.fromRotation(Mat4(), 0.7, Vec3.normalize(Vec3(), Vec3.create(1, 1, 0))),
);

/** each transform of `after` is `move` times the same of `before` */
function expectMoved(before: Float32Array, after: Float32Array, count: number) {
  const expected = Mat4();
  for (let i = 0; i < count; ++i) {
    Mat4.mul(expected, move, Mat4.fromArray(Mat4(), before, i * 16));
    for (let k = 0; k < 16; ++k) expect(after[i * 16 + k]).toBeCloseTo(expected[k], 4);
  }
}

describe("the instances' cell and batch sizes", () => {
  it("are the transform data's, as created and rebuilt, cells above a size of 1, the poses' too", () => {
    const transforms = makeTransforms(30, (i) => Mat4.fromTranslation(Mat4(), Vec3.create(i * 40, 0, 0)));
    const sphere = Sphere3D.create(Vec3(), 1);
    const td = createTransform(transforms, 30, sphere, 100, 8);
    expect([td.instanceCellSize.ref.value, td.instanceBatchSize.ref.value]).toEqual([100, 8]);
    expect(hasInstanceCells(td)).toBe(true);
    // culled as a whole, its levels at full detail: by the size asked for, whatever the grid holds
    updateTransformData(td, sphere, 1, 8);
    expect(hasInstanceCells(td)).toBe(false);
    updateTransformData(td, sphere, 0, 0);
    expect([td.instanceCellSize.ref.value, td.instanceBatchSize.ref.value]).toEqual([0, 0]);
    expect(hasInstanceCells(createIdentityTransform())).toBe(false);
    const positions = new Float32Array(30 * 3).map((_, k) => (k % 3 === 0 ? k * 10 : 0));
    const poses = InstanceTransforms.ofPoses({ positions, count: 30, first: 0, center: Vec3() });
    const posed = createPosesTransform(poses, sphere, 50, 4);
    expect([posed.instanceCellSize.ref.value, posed.instanceBatchSize.ref.value]).toEqual([50, 4]);
    expect(hasInstanceCells(posed)).toBe(true);
  });
});

describe('moveTransformData', () => {
  it('moves the transforms and the grid as a whole, its cells and instance order kept', () => {
    const count = 40;
    const data = createTransform(transforms(count), count, Sphere3D.create(Vec3(), 1), 10, 4);
    const grid = data.instanceGrid.ref.value;
    const transform = data.transform.ref.value.slice();
    const aTransform = data.aTransform.ref.value.slice();
    const aInstance = data.aInstance.ref.value;
    expect(moveTransformData(data, move)).toBe(true);
    expectMoved(transform, data.transform.ref.value, count);
    expectMoved(aTransform, data.aTransform.ref.value, count);
    expect(data.aTransform.ref.value).toBe(data.instanceGrid.ref.value.cellTransform);
    expect(data.instanceGrid.ref.value.layout).toBe(grid);
    expect(data.aInstance.ref.value).toBe(aInstance);
  });

  it('moves the transforms without a grid', () => {
    const count = 5;
    const data = createTransform(transforms(count), count, undefined, 0, 0);
    const transform = data.transform.ref.value.slice();
    const aTransform = data.aTransform.ref.value.slice();
    expect(moveTransformData(data, move)).toBe(true);
    expectMoved(transform, data.transform.ref.value, count);
    expectMoved(aTransform, data.aTransform.ref.value, count);
  });

  it('leaves them where a matrix applies around them', () => {
    const count = 5;
    const data = createTransform(transforms(count), count, Sphere3D.create(Vec3(), 1), 10, 4);
    ValueCell.update(data.matrix, Mat4.fromTranslation(Mat4(), Vec3.create(1, 0, 0)));
    const version = data.aTransform.ref.version;
    expect(moveTransformData(data, move)).toBe(false);
    expect(data.aTransform.ref.version).toBe(version);
  });
});

/** `count` poses of as many particles (scaled, rotated) around `center`, and their matrices (`T(p) R(q) S(s) T(-c)`) */
function makePoses(count: number) {
  const positions = new Float32Array(count * 3);
  const rotations = new Float32Array(count * 4);
  const scales = new Float32Array(count);
  const q = Quat();
  for (let i = 0; i < count; ++i) {
    Vec3.toArray(Vec3.create((i % 7) * 40, Math.floor(i / 7) * 35, (i % 3) * 50), positions, i * 3);
    Quat.setAxisAngle(q, Vec3.normalize(Vec3(), Vec3.create(1, i % 5, 2)), i * 0.3);
    rotations.set(q, i * 4);
    scales[i] = 1 + (i % 4) * 0.5;
  }
  const center = Vec3.create(2, -1, 3);
  const poses: InstancePoses = { kind: 'poses', count, positions, rotations, scales, first: 0, center };
  const matrices = new Float32Array(count * 16);
  const m = Mat4();
  for (let i = 0; i < count; ++i) Mat4.toArray(InstanceTransforms.getMatrix(poses, i, m), matrices, i * 16);
  return { poses, matrices };
}

function rowsOf(t: TransformData) {
  const n = t.instanceCount.ref.value;
  return writeAffineRowsOf(t, n, new Float32Array(n * 12));
}

function expectClose(a: ArrayLike<number>, b: ArrayLike<number>, eps = 1e-3) {
  expect(a.length).toBe(b.length);
  for (let k = 0; k < a.length; ++k) expect(Math.abs(a[k] - b[k])).toBeLessThan(eps);
}

// a particle target's instances as poses read as the same instances as matrices
describe('createPosesTransform', () => {
  const sphere = Sphere3D.create(Vec3.create(2, -1, 3), 5);

  it('draws poses as their matrices: the rows, the grid and its order, no matrix made', () => {
    const { poses, matrices } = makePoses(60);
    const a = createPosesTransform(poses, sphere, 40, 8);
    const b = createTransform(matrices, 60, sphere, 40, 8);
    expect(a.aTransform.ref.value.length).toBe(0);
    expect(a.transform.ref.value.length).toBe(0);
    expect(a.instanceGrid.ref.value.cellTransform.length).toBe(0);
    expect(Array.from(a.aInstance.ref.value)).toEqual(Array.from(b.aInstance.ref.value));
    expectClose(a.instanceGrid.ref.value.cellSpheres, b.instanceGrid.ref.value.cellSpheres, 1e-2);
    expectClose(rowsOf(a), rowsOf(b));
    expectClose(getInstanceMatrices(a), getInstanceMatrices(b));
  });

  it('draws a built grid made pending in the instances own order: the rows as the poses', () => {
    const { poses } = makePoses(60);
    const built = createPosesTransform(poses, sphere, 40, 8);
    // a grid of cells: its instances in the cells' order, not serial
    expect(Array.from(built.aInstance.ref.value)).not.toEqual(Array.from({ length: 60 }, (_, i) => i));
    const pending = createPosesTransform(poses, sphere, 40, 8, built, true);
    expect(pending.instanceGrid.ref.value.pending).toBeDefined();
    // the rows as the poses (no order): in theirs
    expect(Array.from(pending.aInstance.ref.value)).toEqual(Array.from({ length: 60 }, (_, i) => i));
    expectClose(rowsOf(pending), rowsOf(createPosesTransform(poses, sphere, 40, 8, undefined, true)));
  });

  it('applies the representation matrix to the poses, and a move to their matrix without touching an instance', () => {
    const { poses, matrices } = makePoses(30);
    const a = createPosesTransform(poses, sphere, 40, 8);
    const b = createTransform(matrices, 30, sphere, 40, 8);
    const matrix = Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(5, 6, 7)),
      Mat4.fromRotation(Mat4(), 0.4, Vec3.create(0, 1, 0)),
    );
    for (const t of [a, b]) {
      ValueCell.update(t.matrix, Mat4.copy(Mat4(), matrix));
      updateTransformData(t, sphere, 40, 8);
    }
    expectClose(rowsOf(a), rowsOf(b));

    const c = createPosesTransform(poses, sphere, 40, 8);
    const d = createTransform(matrices, 30, sphere, 40, 8);
    const move = Mat4.fromTranslation(Mat4(), Vec3.create(100, 0, -50));
    const grid = c.instanceGrid.ref.value;
    expect(moveTransformData(c, move)).toBe(true);
    expect(moveTransformData(d, move)).toBe(true);
    expect(c.poses.ref.value!.positions).toBe(poses.positions);
    expect(c.instanceGrid.ref.value.cellOffsets).toBe(grid.cellOffsets);
    expectClose(rowsOf(c), rowsOf(d));
    expectClose(c.instanceGrid.ref.value.cellSpheres, d.instanceGrid.ref.value.cellSpheres, 1e-2);
  });

  it('turns the poses into matrices for an extra transform and a mirroring matrix', () => {
    const { poses, matrices } = makePoses(20);
    const a = createPosesTransform(poses, sphere, 40, 8);
    const b = createTransform(matrices, 20, sphere, 40, 8);
    const mirror = Mat4.fromScaling(Mat4(), Vec3.create(-1, 1, 1));
    for (const t of [a, b]) {
      ValueCell.update(t.matrix, Mat4.copy(Mat4(), mirror));
      updateTransformData(t, sphere, 40, 8);
    }
    expect(a.poses.ref.value).toBeUndefined();
    expect(a.hasReflection.ref.value).toBe(b.hasReflection.ref.value);
    expectClose(rowsOf(a), rowsOf(b));
  });
});
