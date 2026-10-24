/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { DragFrame } from '@molstar/graphics/util/gizmo/gestures';
import { describe, expect, it } from 'vitest';
import { type CutObject, coneAngle, coneScale } from '../../cuts.js';
import { applyToCut } from '../binding.js';
import { ConeAngleHandle } from '../geometry.js';

// What a drag on a cut's handle makes of the cut: the gizmo library's gestures on its pose, back to a clip object
// (only what changed), and the cone's angle (the binding's own).

const frame: DragFrame = { right: Vec3.create(1, 0, 0), up: Vec3.create(0, 1, 0), perPixel: 0.5 };
const close = (a: ArrayLike<number>, b: readonly number[]) => b.forEach((x, i) => expect(a[i]).toBeCloseTo(x, 6));

const cut = (
  type: CutObject['type'],
  scale = Vec3.create(10, 10, 10),
  rotation = { axis: Vec3.create(1, 0, 0), angle: 0 },
): CutObject => ({
  type,
  invert: false,
  position: Vec3.create(0, 0, 0),
  rotation,
  scale,
  transform: undefined as never,
});

describe('the clipping gizmos binding', () => {
  it('moves a cut and leaves its rotation as it was (axis and angle)', () => {
    const turned = cut('cube', undefined, { axis: Vec3.create(0, 0, 3), angle: -30 });
    const next = applyToCut(turned, { kind: 'move' }, frame, 20, 0);
    close(next.position, [10, 0, 0]);
    expect(next.rotation).toBe(turned.rotation);
    expect(next.scale).toBe(turned.scale);
  });

  it('turns a cut whose stored axis is not unit length', () => {
    const from = cut('cube', undefined, { axis: Vec3.create(0, 0, 2), angle: 30 });
    const next = applyToCut(from, { kind: 'turn' }, frame, 0, 0);
    expect(next.rotation.angle).toBeCloseTo(30, 4);
    close(next.rotation.axis, [0, 0, 1]);
  });

  it("drags a turned box's face along its own axis", () => {
    const turned = cut('cube', undefined, { axis: Vec3.create(0, 0, 1), angle: 90 });
    // its local +X is world +Y: dragging up the screen grows it
    const next = applyToCut(turned, { kind: 'faces', faces: [{ axis: 0, sign: 1 }] }, frame, 0, -8);
    expect(next.scale[0]).toBeCloseTo(14, 6);
    close(next.position, [0, 2, 0]);
  });

  it('widens a cylinder by its rim: x and z, its height kept', () => {
    const cylinder = cut('cylinder', Vec3.create(10, 20, 10));
    const next = applyToCut(cylinder, { kind: 'scale', axes: [0, 2] }, { ...frame, grab: [20, 0] }, 10, 0);
    close(next.scale, [15, 20, 15]);
  });

  it("widens a cone: its base rim's radius goes as the angle's tangent", () => {
    const cone = cut('infiniteCone', coneScale(30));
    const next = applyToCut(cone, { kind: 'custom', id: ConeAngleHandle }, { ...frame, grab: [20, 0] }, 20, 0);
    expect(Math.tan((coneAngle(next.scale) * Math.PI) / 180)).toBeCloseTo(2 * Math.tan(Math.PI / 6), 6);
  });

  it('does not size a plane', () => {
    const plane = cut('plane');
    expect(applyToCut(plane, { kind: 'scale' }, { ...frame, grab: [20, 0] }, 20, 0)).toBe(plane);
  });
});
