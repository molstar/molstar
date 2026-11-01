/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { Sphere3D } from '@molstar/core/math/geometry/primitives/sphere3d';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Cylinders, CylinderStride } from '../cylinders/cylinders.js';
import { Lines, LineStride } from '../lines/lines.js';

/** `count` items of `stride` floats, each with a start at 0 and an end at 4, on a line of growing distance */
function createItems(count: number, stride: number) {
  const data = new Float32Array(count * stride);
  for (let i = 0; i < count; ++i) {
    const o = i * stride;
    data.set([i, Math.sin(i) * 5, Math.cos(i) * 5, i, i + 1, Math.sin(i + 1) * 5, Math.cos(i + 1) * 5], o);
  }
  return data;
}

function expectEncloses(sphere: Sphere3D, data: Float32Array, count: number, stride: number) {
  for (let i = 0; i < count; ++i) {
    for (const offset of [0, 4]) {
      const p = Vec3.fromArray(Vec3(), data, i * stride + offset);
      expect(Vec3.distance(p, sphere.center)).toBeLessThanOrEqual(sphere.radius * (1 + 1e-5));
    }
  }
}

// the item counts choose the helper and the exact extrema, not the vertices items once had
describe('bounding spheres of geometries with data per item', () => {
  for (const count of [1, 16, 17, 49, 98, 99, 2000, 20_000]) {
    it(`encloses all ${count} cylinders and lines`, () => {
      const cylinders = createItems(count, CylinderStride);
      expectEncloses(Cylinders.create(cylinders, count).boundingSphere, cylinders, count, CylinderStride);
      const lines = createItems(count, LineStride);
      expectEncloses(Lines.create(lines, count).boundingSphere, lines, count, LineStride);
    });
  }

  it('keeps few cylinders as exact extrema', () => {
    const data = createItems(40, CylinderStride);
    // two spheres of 40 positions each, merged: their points, not directions' samples
    const { extrema } = Cylinders.create(data, 40).boundingSphere;
    for (const p of extrema ?? []) {
      expect(Number.isFinite(p[0])).toBe(true);
    }
    expect(extrema?.length ?? 0).toBeGreaterThan(0);
  });
});
