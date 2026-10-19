/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Mat4, Vec3 } from '../../linear-algebra.js';
import { Box3D } from '../primitives/box3d.js';

describe('Box3D', () => {
  it('fromPositionArray', () => {
    const positions = Float32Array.of(1, 2, 3, -4, 5, 0.5, 2, -1, 9, 100, 100, 100);
    const box = Box3D.fromPositionArray(Box3D(), positions, 3);
    expect(Array.from(box.min)).toEqual([-4, -1, 0.5]);
    expect(Array.from(box.max)).toEqual([2, 5, 9]);
    expect(Box3D.isEmpty(Box3D.fromPositionArray(Box3D(), positions, 0))).toBe(true);
  });

  it('intersect and isEmpty', () => {
    const a = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(2, 2, 2));
    const b = Box3D.create(Vec3.create(1, -1, 1), Vec3.create(3, 1, 5));
    const ab = Box3D.intersect(Box3D(), a, b);
    expect(Array.from(ab.min)).toEqual([1, 0, 1]);
    expect(Array.from(ab.max)).toEqual([2, 1, 2]);
    expect(Box3D.isEmpty(ab)).toBe(false);
    const far = Box3D.create(Vec3.create(5, 5, 5), Vec3.create(6, 6, 6));
    expect(Box3D.isEmpty(Box3D.intersect(Box3D(), a, far))).toBe(true);
    expect(Box3D.isEmpty(Box3D.setEmpty(Box3D()))).toBe(true);
    // a point is not empty
    expect(Box3D.isEmpty(Box3D.create(Vec3.create(1, 1, 1), Vec3.create(1, 1, 1)))).toBe(false);
  });

  it('transform', () => {
    const box = Box3D.create(Vec3.create(0, 0, 0), Vec3.create(1, 2, 3));
    const m = Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.create(0, 0, 1));
    Mat4.setTranslation(m, Vec3.create(10, 0, 0));
    const t = Box3D.transform(Box3D(), box, m);
    expect(Array.from(t.min).map((x) => Math.round(x * 1e6) / 1e6)).toEqual([8, 0, 0]);
    expect(Array.from(t.max).map((x) => Math.round(x * 1e6) / 1e6)).toEqual([10, 1, 3]);
  });
});
