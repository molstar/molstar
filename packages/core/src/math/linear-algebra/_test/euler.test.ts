/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { describe, expect, it } from 'vitest';
import { Mat4 } from '../3d/mat4.js';
import { Euler } from '../3d/euler.js';
import { Quat } from '../3d/quat.js';

const t = [
  [Euler.create(0, 0, 0), 'XYZ'],
  [Euler.create(1, 0, 0), 'XYZ'],
  [Euler.create(0, 1, 0), 'ZYX'],
] as const;

describe('Euler', () => {
  it('fromMat4', () => {
    for (const [e, o] of t) {
      const m = Mat4.fromEuler(Mat4(), e, o);
      const e2 = Euler.fromMat4(Euler(), m, o);
      const m2 = Mat4.fromEuler(Mat4(), e2, o);
      expect(Mat4.areEqual(m, m2, 0.0001)).toBe(true);
    }
  });

  it('set and fromVec3 write each angle', () => {
    expect(Array.from(Euler.set(Euler(), 1, 2, 3))).toEqual([1, 2, 3]);
    expect(Array.from(Euler.fromVec3(Euler(), [4, 5, 6] as never))).toEqual([4, 5, 6]);
  });

  it('fromQuat', () => {
    for (const [e, o] of t) {
      const q = Quat.fromEuler(Quat(), e, o);
      const e2 = Euler.fromQuat(Euler(), q, o);
      const q2 = Quat.fromEuler(Quat(), e2, o);
      expect(Quat.equals(q, q2)).toBe(true);
    }
  });
});
