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
import { createIdentityTransform, createTransform, moveTransformData, updateTransformData } from '../transform-data.js';

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
