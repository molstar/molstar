/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { OrderedSet } from '../../../data/int/ordered-set.js';
import { Vec3 } from '../../linear-algebra/3d/vec3.js';
import { getBoundary } from '../boundary.js';
import { BoundaryHelper } from '../boundary-helper.js';
import { Sphere3D } from '../primitives/sphere3d.js';

/** `n` points on a jittered line along x (elongated: tight spheres need the points), radii 1 to 2 */
function points(n: number, offset = 0) {
  const x: number[] = [];
  const y: number[] = [];
  const z: number[] = [];
  const radius: number[] = [];
  let seed = 7;
  const rand = () => {
    seed = (seed * 16807) % 2147483647;
    return seed / 2147483647;
  };
  for (let i = 0; i < n; ++i) {
    x.push(offset + i * 1.5);
    y.push(rand() - 0.5);
    z.push(rand() - 0.5);
    radius.push(1 + rand());
  }
  return { x, y, z, radius, indices: OrderedSet.ofBounds(0, n) };
}

/** Whether `s` covers every ball of `data` (within float tolerance) */
function coversBalls(s: Sphere3D, data: ReturnType<typeof points>, withRadius = true) {
  for (let i = 0; i < data.x.length; ++i) {
    const d = Vec3.distance(s.center, Vec3.create(data.x[i], data.y[i], data.z[i]));
    if (d + (withRadius ? data.radius[i] : 0) > s.radius + 1e-4) return false;
  }
  return true;
}

describe('boundary extrema', () => {
  it('keeps few points with radii as extrema, with their largest radius', () => {
    const data = points(6);
    const { sphere } = getBoundary(data);
    expect(sphere.extrema?.length).toBe(6);
    expect(sphere.extremaRadius).toBe(Math.max(...data.radius));
    expect(coversBalls(sphere, data)).toBe(true);
  });

  it('keeps few points without radii as extrema, without a radius', () => {
    const data = points(6);
    const { sphere } = getBoundary({ ...data, radius: undefined });
    expect(sphere.extrema?.length).toBe(6);
    expect(sphere.extremaRadius).toBeUndefined();
  });

  it('keeps no more extrema than distinct points', () => {
    // 300 points, many of them extreme in several directions
    const data = points(300);
    const { sphere } = getBoundary({ ...data, radius: undefined });
    const extrema = sphere.extrema!;
    expect(extrema.length).toBeLessThanOrEqual(98);
    const keys = new Set(extrema.map((e) => `${e[0]},${e[1]},${e[2]}`));
    expect(keys.size).toBe(extrema.length);
  });

  it('a union of spheres of few points with radii covers their balls', () => {
    const a = points(6);
    const b = points(5, 30);
    const sa = getBoundary(a).sphere;
    const sb = getBoundary(b).sphere;
    const helper = new BoundaryHelper('98');
    helper.reset();
    helper.includeSphere(sa);
    helper.includeSphere(sb);
    helper.finishedIncludeStep();
    helper.radiusSphere(sa);
    helper.radiusSphere(sb);
    const union = helper.getSphere();
    expect(coversBalls(union, a)).toBe(true);
    expect(coversBalls(union, b)).toBe(true);
    // tighter than the union of the two spheres
    const loose = Sphere3D.expandBySphere(Sphere3D(), sa, sb);
    expect(union.radius).toBeLessThan(loose.radius);
  });
});
