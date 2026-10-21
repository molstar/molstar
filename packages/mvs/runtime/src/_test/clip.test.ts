/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { Plane3D } from '@molstar/core/math/geometry/primitives/plane3d';
import { Clip } from '@molstar/core/util/clip';
import { clippingForNode } from '../load-helpers.js';

function planeClip(normal: [number, number, number]) {
  const node = {
    kind: 'representation',
    params: { type: 'cartoon' },
    children: [{ kind: 'clip', params: { type: 'plane', point: [0, 0, 0], normal, variant: 'pixel', invert: false } }],
  } as any;
  const props = clippingForNode(node)!;
  return { props, plane: Clip.getPlane(Plane3D(), Clip.getClip(props).objects, 0) };
}

/** Positive side of the plane is the clipped half-space (`Clip.getPlane`). */
function isClipped(plane: Plane3D, p: Vec3) {
  return Vec3.dot(plane.normal, p) + plane.constant > 0;
}

describe('clippingForNode, plane', () => {
  it('normal down turns up by 180°', () => {
    const { props, plane } = planeClip([0, -1, 0]);
    expect(props.objects[0].rotation.angle).toBeCloseTo(180);
    expect(Vec3.magnitude(props.objects[0].rotation.axis)).toBeCloseTo(1);
    expect(plane.normal[1]).toBeCloseTo(-1);
    expect(isClipped(plane, Vec3.create(0, -1, 0))).toBe(true);
    expect(isClipped(plane, Vec3.create(0, 1, 0))).toBe(false);
  });

  it('normal up keeps 0°', () => {
    const { props, plane } = planeClip([0, 2, 0]);
    expect(props.objects[0].rotation.angle).toBeCloseTo(0);
    expect(plane.normal[1]).toBeCloseTo(1);
    expect(isClipped(plane, Vec3.create(0, 1, 0))).toBe(true);
    expect(isClipped(plane, Vec3.create(0, -1, 0))).toBe(false);
  });

  it('normal +z turns by 90°', () => {
    const { props, plane } = planeClip([0, 0, 1]);
    expect(props.objects[0].rotation.angle).toBeCloseTo(90);
    expect(plane.normal[2]).toBeCloseTo(1);
    expect(isClipped(plane, Vec3.create(0, 0, 1))).toBe(true);
    expect(isClipped(plane, Vec3.create(0, 0, -1))).toBe(false);
  });
});
