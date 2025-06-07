/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4 } from '../3d/mat4';
import { Vec3 } from '../3d/vec3';
import { describe, it, expect } from 'vitest';

describe('Mat4', () => {
    it('permutation', () => {
        expect(Mat4.areEqual(Mat4.fromPermutation(Mat4(), [0, 1, 2, 3]), Mat4.identity(), 1e-6)).toBe(true);

        expect(Mat4.areEqual(Mat4.fromPermutation(Mat4(), [1, 0, 2, 3]), Mat4.ofRows([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]), 1e-6)).toBe(true);

        const perm = Mat4.fromPermutation(Mat4(), [1, 2, 0, 3]);

        expect(Vec3.transformMat4(Vec3(), Vec3.create(1, 2, 3), perm)).toEqual(Vec3.create(2, 3, 1));
    });
});