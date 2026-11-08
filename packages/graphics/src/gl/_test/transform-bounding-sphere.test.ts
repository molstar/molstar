/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { calculateTransformBoundingSphere } from '../renderable/util.js';

/** `count` rigid transforms, rotated and spread over a slab */
function transforms(count: number) {
  const array = new Float32Array(count * 16);
  const m = Mat4();
  const q = Quat();
  let seed = 1;
  const random = () => ((seed = (seed * 16807) % 2147483647) - 1) / 2147483646;
  for (let i = 0; i < count; ++i) {
    Quat.normalize(q, Quat.create(random() - 0.5, random() - 0.5, random() - 0.5, random() - 0.5));
    Mat4.fromQuat(m, q);
    Mat4.setTranslation(m, Vec3.create(random() * 1000, random() * 400, random() * 50));
    Mat4.toArray(m, array, i * 16);
  }
  return array;
}

describe('calculateTransformBoundingSphere', () => {
  const invariant = Sphere3D.create(Vec3.create(3, -2, 5), 4);
  const count = 2000;
  const array = transforms(count);

  it('encloses every instance, many of them (the box sphere)', () => {
    const sphere = calculateTransformBoundingSphere(invariant, array, count, 0);
    const m = Mat4();
    const s = Sphere3D();
    for (let i = 0; i < count; ++i) {
      Sphere3D.transform(s, invariant, Mat4.fromArray(m, array, i * 16));
      expect(Vec3.distance(s.center, sphere.center) + s.radius).toBeLessThanOrEqual(sphere.radius * (1 + 1e-6));
    }
  });

  it('is about as tight from 1000 instances (the box sphere) as below (the helper)', () => {
    const box = calculateTransformBoundingSphere(invariant, array, 1000, 0);
    const helper = calculateTransformBoundingSphere(invariant, array, 999, 0);
    expect(box.radius).toBeGreaterThan(helper.radius * 0.95);
    expect(box.radius).toBeLessThan(helper.radius * 1.05);
  });

  it('keeps few instances spread wide tight', () => {
    // a dozen spread over a shell of radius 1500: the box's center would be off
    const shell = new Float32Array(12 * 16);
    const m = Mat4();
    for (let i = 0; i < 12; ++i) {
      const t = (i * 2.4) % (2 * Math.PI);
      const z = 1 - (2 * i + 1) / 12;
      const r = Math.sqrt(1 - z * z);
      Mat4.fromTranslation(m, Vec3.create(1500 * r * Math.cos(t), 1500 * r * Math.sin(t), 1500 * z));
      Mat4.toArray(m, shell, i * 16);
    }
    const sphere = calculateTransformBoundingSphere(invariant, shell, 12, 0);
    expect(sphere.radius).toBeLessThan((1500 + Vec3.magnitude(invariant.center) + invariant.radius) * 1.02);
  });

  it('keeps a single instance exact', () => {
    const one = calculateTransformBoundingSphere(invariant, array, 1, 0);
    const expected = Sphere3D.transform(Sphere3D(), invariant, Mat4.fromArray(Mat4(), array, 0));
    expect(Vec3.distance(one.center, expected.center)).toBeLessThan(1e-4);
    expect(one.radius).toBeCloseTo(expected.radius, 4);
  });
});
