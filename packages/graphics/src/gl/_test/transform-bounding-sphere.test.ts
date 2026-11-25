/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { InstanceTransforms } from '@molstar/core/math/geometry/instance-transforms';
import { ValueCell } from '@molstar/core/util/value-cell';
import {
  calculateTransformBoundingSphere,
  calculateTransformDataBoundingSphere,
  setInstanceSphereFast,
  setTransformDataCenters,
} from '../renderable/util.js';
import { createTransform } from '../../geo/geometry/transform-data.js';

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

describe('calculateTransformDataBoundingSphere', () => {
  const count = 2000;
  const values = () => ({
    aTransform: ValueCell.create(transforms(count)),
    instanceCount: ValueCell.create(count),
  });

  it('is the instances sphere, a radius changed alone kept from the pass before', () => {
    const data = values();
    const small = Sphere3D.create(Vec3.create(3, -2, 5), 4);
    const large = Sphere3D.create(Vec3.create(3, -2, 5), 9);
    const a = calculateTransformDataBoundingSphere(small, data);
    expect(a).toEqual(calculateTransformBoundingSphere(small, data.aTransform.ref.value, count, 0));
    // the transforms unchanged: the same center, the distance plus the radius
    const writes = vi.spyOn(InstanceTransforms, 'writePoints');
    writes.mockClear();
    const b = calculateTransformDataBoundingSphere(large, data);
    expect(writes).not.toHaveBeenCalled();
    expect(b).toEqual(calculateTransformBoundingSphere(large, data.aTransform.ref.value, count, 0));
    expect(b.radius - a.radius).toBeCloseTo(5, 4);
    writes.mockRestore();
  });

  it('passes over the instances again where they changed', () => {
    const data = values();
    const invariant = Sphere3D.create(Vec3.create(3, -2, 5), 4);
    const before = calculateTransformDataBoundingSphere(invariant, data);
    const array = data.aTransform.ref.value;
    for (let i = 0; i < count; ++i) array[i * 16 + 12] += 500;
    ValueCell.update(data.aTransform, array);
    const after = calculateTransformDataBoundingSphere(invariant, data);
    expect(after.center[0] - before.center[0]).toBeCloseTo(500, 2);
    expect(after).toEqual(calculateTransformBoundingSphere(invariant, array, count, 0));
  });

  it("keeps a few instances' sphere (the helper's) for the same transforms and invariant sphere", () => {
    const few = 40;
    const data = { aTransform: ValueCell.create(transforms(few)), instanceCount: ValueCell.create(few) };
    const invariant = Sphere3D.create(Vec3.create(3, -2, 5), 4);
    const points = vi.spyOn(InstanceTransforms, 'getPoint');
    const a = calculateTransformDataBoundingSphere(invariant, data);
    expect(a).toEqual(calculateTransformBoundingSphere(invariant, data.aTransform.ref.value, few, 0));
    points.mockClear();
    // the same: no pass; another radius (not additive for the helper's sphere) or other transforms: one
    expect(calculateTransformDataBoundingSphere(Sphere3D.clone(invariant), data)).toEqual(a);
    expect(points).not.toHaveBeenCalled();
    const larger = Sphere3D.create(Vec3.create(3, -2, 5), 9);
    expect(calculateTransformDataBoundingSphere(larger, data)).toEqual(
      calculateTransformBoundingSphere(larger, data.aTransform.ref.value, few, 0),
    );
    ValueCell.update(data.aTransform, data.aTransform.ref.value);
    points.mockClear();
    calculateTransformDataBoundingSphere(larger, data);
    expect(points).toHaveBeenCalled();
    points.mockRestore();
  });
});

describe('the instances bounding sphere asked fast', () => {
  it("is a few instances' box sphere, seeded without a pass", () => {
    const few = 40;
    const data = { aTransform: ValueCell.create(transforms(few)), instanceCount: ValueCell.create(few) };
    setInstanceSphereFast(data, true);
    const invariant = Sphere3D.create(Vec3.create(3, -2, 5), 4);
    const boxSphere = calculateTransformBoundingSphere(invariant, data.aTransform.ref.value, few, 0);
    const sphere = calculateTransformDataBoundingSphere(invariant, data);
    // a box sphere holds the helper's sphere of the same instances: no smaller
    expect(sphere.radius).toBeGreaterThanOrEqual(boxSphere.radius - 1e-6);
    const points = vi.spyOn(InstanceTransforms, 'getPoint');
    const writes = vi.spyOn(InstanceTransforms, 'writePoints');
    ValueCell.update(data.aTransform, data.aTransform.ref.value);
    setTransformDataCenters(data, { invariantCenter: invariant.center, center: Vec3.create(1, 2, 3), distance: 10 });
    const seeded = calculateTransformDataBoundingSphere(invariant, data);
    expect(points).not.toHaveBeenCalled();
    expect(writes).not.toHaveBeenCalled();
    expect(Array.from(seeded.center)).toEqual([1, 2, 3]);
    expect(seeded.radius).toBe(14);
    points.mockRestore();
    writes.mockRestore();
  });
});

describe('the instances bounding sphere from their grid', () => {
  it('needs no pass over the instances after the grid built from them', () => {
    const count = 2000;
    const array = transforms(count);
    const invariant = Sphere3D.create(Vec3.create(3, -2, 5), 4);
    const data = createTransform(array.slice(), count, invariant, 50, 64);
    // the padded sphere (the same center) as the geometry's values have it
    const padded = Sphere3D.create(Vec3.create(3, -2, 5), 6);
    const writes = vi.spyOn(InstanceTransforms, 'writePoints');
    writes.mockClear();
    const sphere = calculateTransformDataBoundingSphere(padded, data);
    expect(writes).not.toHaveBeenCalled();
    const expected = calculateTransformBoundingSphere(padded, array, count, 0);
    expect(Vec3.distance(sphere.center, expected.center)).toBeLessThan(1e-4);
    expect(sphere.radius).toBeCloseTo(expected.radius, 3);
    writes.mockRestore();
  });
});
