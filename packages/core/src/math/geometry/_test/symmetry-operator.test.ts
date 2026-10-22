/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { Mat4 } from '../../linear-algebra/3d/mat4.js';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { SymmetryOperator } from '../symmetry-operator.js';

const op = (m: Mat4) => SymmetryOperator.create('op', m);
const translation = (x: number, y: number, z: number) => Mat4.fromTranslation(Mat4(), Vec3.create(x, y, z));
const motion = () =>
  Mat4.setTranslation(Mat4.fromRotation(Mat4(), Math.PI / 3, Vec3.create(1, 2, 3)), Vec3.create(4, 5, 6));
const moved = (m: Mat4, a: SymmetryOperator) => op(Mat4.mul(Mat4(), m, a.matrix));

describe('SymmetryOperator.rigidMotionTest', () => {
  const ops = [op(Mat4.identity()), op(translation(10, 0, 0)), op(Mat4.fromRotation(Mat4(), 1, Vec3.create(0, 1, 0)))];

  it('accepts the identity: equal operators, not moved', () => {
    const test = SymmetryOperator.rigidMotionTest();
    for (const a of ops) expect(test(a, a)).toBe(true);
    expect(test.moved).toBe(false);
    // after an identity first pair, a moved one is another motion
    expect(test(ops[1], moved(motion(), ops[1]))).toBe(false);
  });

  it('accepts every operator after one common rigid motion', () => {
    const m = motion();
    const test = SymmetryOperator.rigidMotionTest();
    for (const a of ops) expect(test(a, moved(m, a))).toBe(true);
    expect(test.moved).toBe(true);
  });

  it('rejects a second motion, an unmoved operator after a moved one, and a non-rigid motion', () => {
    const test = SymmetryOperator.rigidMotionTest();
    expect(test(ops[0], moved(motion(), ops[0]))).toBe(true);
    expect(test(ops[1], moved(translation(0, 1, 0), ops[1]))).toBe(false);
    expect(test(ops[2], ops[2])).toBe(false);

    // a uniform scaling keeps no distance (`create` warns about it)
    vi.spyOn(console, 'warn').mockImplementation(() => {});
    const scaling = SymmetryOperator.rigidMotionTest();
    expect(scaling(ops[1], moved(Mat4.fromScaling(Mat4(), Vec3.create(2, 2, 2)), ops[1]))).toBe(false);
  });
});
