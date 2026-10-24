/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { Box3D } from '../primitives/box3d.js';
import { Sphere3D } from '../primitives/sphere3d.js';

describe('Sphere3D', () => {
  const withExtrema = () => {
    const s = Sphere3D.create(Vec3.create(0, 0, 0), 1);
    Sphere3D.setExtrema(s, [Vec3.create(-1, 0, 0), Vec3.create(1, 0, 0)]);
    return s;
  };

  it('copy from a sphere without extrema clears them (perf-audit-2 02/B3)', () => {
    const out = withExtrema();
    Sphere3D.copy(out, Sphere3D.create(Vec3.create(10, 0, 0), 2));
    expect(out.extrema).toBeUndefined();
    expect(Sphere3D.hasExtrema(out)).toBe(false);
  });

  it('expandBySphere by a sphere without extrema clears them', () => {
    const out = withExtrema();
    Sphere3D.expandBySphere(out, out, Sphere3D.create(Vec3.create(10, 0, 0), 2));
    expect(out.radius).toBe(12);
    expect(out.extrema).toBeUndefined();
    // a box from stale extrema would bound the old points only
    const box = Box3D.fromSphere3D(Box3D(), out);
    expect(box.max[0]).toBeGreaterThanOrEqual(12);
  });

  it('copy and expandBySphere keep extrema when both have them', () => {
    const out = Sphere3D.copy(Sphere3D(), withExtrema());
    expect(out.extrema?.length).toBe(2);
    const by = Sphere3D.create(Vec3.create(3, 0, 0), 1);
    Sphere3D.setExtrema(by, [Vec3.create(4, 0, 0)]);
    Sphere3D.expandBySphere(out, out, by);
    expect(out.extrema?.length).toBe(3);
  });
});
