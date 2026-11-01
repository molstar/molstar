/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
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

  it('expand keeps the extrema and grows the balls around them', () => {
    const s = withExtrema();
    const out = Sphere3D.expand(Sphere3D(), s, 0.5);
    expect(out.radius).toBe(1.5);
    expect(out.extrema).toEqual(s.extrema);
    expect(out.extremaRadius).toBe(0.5);
    // in place too
    Sphere3D.expand(out, out, 0.25);
    expect(out.extrema?.length).toBe(2);
    expect(out.extremaRadius).toBe(0.75);
  });

  it('transform and scale scale the radius around the extrema', () => {
    const s = Sphere3D.setExtrema(Sphere3D.create(Vec3(), 2), [Vec3.create(-1, 0, 0), Vec3.create(1, 0, 0)], 1);
    expect(Sphere3D.scale(Sphere3D(), s, 2).extremaRadius).toBe(2);
    expect(Sphere3D.translate(Sphere3D(), s, Vec3.create(1, 0, 0)).extremaRadius).toBe(1);
    expect(Sphere3D.clone(s).extremaRadius).toBe(1);
  });

  it('setExtrema leaves no radius from an earlier set', () => {
    const s = Sphere3D.setExtrema(Sphere3D.create(Vec3(), 2), [Vec3.create(-1, 0, 0), Vec3.create(1, 0, 0)], 1);
    Sphere3D.setExtrema(s, [Vec3.create(-2, 0, 0), Vec3.create(2, 0, 0)]);
    expect(s.extremaRadius).toBeUndefined();
  });

  it('includes and boxes take the balls around the extrema', () => {
    const s = Sphere3D.setExtrema(Sphere3D.create(Vec3(), 2), [Vec3.create(-1, 0, 0), Vec3.create(1, 0, 0)], 1);
    expect(Sphere3D.includes(Sphere3D.create(Vec3(), 1.5), s)).toBe(false);
    expect(Sphere3D.includes(Sphere3D.create(Vec3(), 2), s)).toBe(true);
    // 14 extrema or more are boxed: the box grown by the radius
    const many = Array.from({ length: 14 }, (_, i) => Vec3.create(i % 2 ? 1 : -1, 0, 0));
    Sphere3D.setExtrema(s, many, 1);
    const box = Box3D.fromSphere3D(Box3D(), s);
    expect(Array.from(box.min)).toEqual([-2, -1, -1]);
    expect(Array.from(box.max)).toEqual([2, 1, 1]);
  });
});
