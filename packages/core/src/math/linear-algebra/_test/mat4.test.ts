/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4 } from '../3d/mat4.js';
import { Vec3 } from '../3d/vec3.js';

describe('Mat4', () => {
  it('permutation', () => {
    expect(Mat4.areEqual(Mat4.fromPermutation(Mat4(), [0, 1, 2, 3]), Mat4.identity(), 1e-6)).toBe(true);

    expect(
      Mat4.areEqual(
        Mat4.fromPermutation(Mat4(), [1, 0, 2, 3]),
        Mat4.ofRows([
          [0, 1, 0, 0],
          [1, 0, 0, 0],
          [0, 0, 1, 0],
          [0, 0, 0, 1],
        ]),
        1e-6,
      ),
    ).toBe(true);

    const perm = Mat4.fromPermutation(Mat4(), [1, 2, 0, 3]);

    expect(Vec3.transformMat4(Vec3(), Vec3.create(1, 2, 3), perm)).toEqual(Vec3.create(2, 3, 1));
  });

  it('isRigid', () => {
    const rotation = Mat4.fromRotation(Mat4(), 0.7, Vec3.normalize(Vec3(), Vec3.create(1, 2, 3)));
    const rigid = Mat4.mul(Mat4(), Mat4.fromTranslation(Mat4(), Vec3.create(10, -5, 3)), rotation);
    expect(Mat4.isRigid(Mat4.identity())).toBe(true);
    expect(Mat4.isRigid(rigid)).toBe(true);
    // a reflection keeps distances
    expect(Mat4.isRigid(Mat4.fromScaling(Mat4(), Vec3.create(-1, 1, 1)))).toBe(true);
    // scaling, and a shear that keeps the volume (isRotationAndTranslation accepts it), do not
    expect(Mat4.isRigid(Mat4.fromScaling(Mat4(), Vec3.create(1.01, 1, 1)))).toBe(false);
    const shear = Mat4.identity();
    shear[4] = 0.5;
    expect(Mat4.isRotationAndTranslation(shear)).toBe(true);
    expect(Mat4.isRigid(shear)).toBe(false);
    // projective
    const projective = Mat4.clone(rigid);
    projective[3] = 0.1;
    expect(Mat4.isRigid(projective)).toBe(false);
    // tolerance
    const nearly = Mat4.fromScaling(Mat4(), Vec3.create(1 + 1e-5, 1, 1));
    expect(Mat4.isRigid(nearly)).toBe(false);
    expect(Mat4.isRigid(nearly, 1e-4)).toBe(true);
  });
});
