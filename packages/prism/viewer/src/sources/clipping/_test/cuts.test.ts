/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import {
  AroundRoom,
  reshape,
  coneAngle,
  coneScale,
  LocalUp,
  MinPlacedRadius,
  placeCut,
  rotated,
  rotationTo,
} from '../cuts.js';
import { applyToCut } from '../gizmo/binding.js';

// The cut geometry the area and the gizmos share: rotations, a cone's angle, where a new cut goes.

const close = (a: ArrayLike<number>, b: readonly number[], digits = 6) =>
  b.forEach((x, i) => expect(a[i]).toBeCloseTo(x, digits));

describe('cuts', () => {
  it('turns +Y onto any direction, opposite included', () => {
    for (const to of [Vec3.create(1, 0, 0), Vec3.create(0, 0, -1), Vec3.create(0, -1, 0), Vec3.create(0, 1, 0)]) {
      close(rotated(rotationTo(LocalUp, to), LocalUp), [...to]);
    }
  });

  it("keeps a cone's half-angle through its scale", () => {
    for (const a of [1, 25, 45, 89]) expect(coneAngle(coneScale(a))).toBeCloseTo(a, 6);
  });

  it('puts a closed shape around a small target, with room', () => {
    const drawn = Sphere3D.create(Vec3.create(0, 0, 0), 30);
    const ligand = Sphere3D.create(Vec3.create(5, 5, 5), 10);
    const sphere = placeCut('sphere', ligand, Vec3.create(0, 0, 1), drawn);
    close(sphere.position, [5, 5, 5]);
    expect(sphere.scale[0]).toBeCloseTo(2 * (10 + AroundRoom), 6);
    // never smaller than a cut one can see and grab
    const atom = placeCut('cube', Sphere3D.create(Vec3.create(5, 5, 5), 1), Vec3.create(0, 0, 1), drawn);
    expect(atom.scale[0]).toBeCloseTo(2 * MinPlacedRadius, 6);
  });

  it('opens a window into a large target, toward the viewer', () => {
    const polymer = Sphere3D.create(Vec3.create(0, 0, 0), 25);
    const sphere = placeCut('sphere', polymer, Vec3.create(0, 0, 1), polymer);
    close(sphere.position, [0, 0, 0.85 * 25]);
    // a plane through the center, facing the viewer
    const plane = placeCut('plane', polymer, Vec3.create(0, 0, 1), polymer);
    close(plane.position, [0, 0, 0]);
    close(rotated(plane.rotation, LocalUp), [0, 0, 1]);
  });

  it('turns a cut whose stored axis is not unit length', () => {
    const from = {
      type: 'cube' as const,
      invert: false,
      position: Vec3(),
      rotation: { axis: Vec3.create(0, 0, 2), angle: 30 },
      scale: Vec3.create(1, 1, 1),
      transform: undefined as never,
    };
    const frame = { right: Vec3.create(1, 0, 0), up: Vec3.create(0, 1, 0), perPixel: 1 };
    const next = applyToCut(from, { kind: 'turn' }, frame, 0, 0);
    expect(next.rotation.angle).toBeCloseTo(30, 4);
  });

  it('reshapes a cut, its place, direction, size and flip carried over', () => {
    const plane = placeCut('plane', Sphere3D.create(Vec3.create(1, 2, 3), 20), Vec3.create(1, 0, 0));
    const flipped = { ...plane, invert: true };
    // a plane has no size: the new one takes the given one; its normal becomes the cylinder's axis
    const cylinder = reshape(flipped, 'cylinder', 16);
    close(cylinder.position, [1, 2, 3]);
    close(rotated(cylinder.rotation, LocalUp), [1, 0, 0]);
    expect(cylinder.scale[0]).toBeCloseTo(16, 6);
    expect(cylinder.invert).toBe(true);
    // and the cone's opening
    const cone = reshape(cylinder, 'infiniteCone', 16);
    close(rotated(cone.rotation, Vec3.create(0, 0, -1)), [1, 0, 0]);
    // a sphere's diameter is a box's edge
    const box = reshape({ ...cylinder, type: 'sphere', scale: Vec3.create(9, 9, 9) }, 'cube', 16);
    close(box.scale, [9, 9, 9]);
    expect(reshape(box, 'cube', 16)).toBe(box);
  });
});
