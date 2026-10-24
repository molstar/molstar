/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import {
  type DragFrame,
  type GizmoPose,
  MaxSlideGain,
  MinExtent,
  inPlane,
  move,
  poseGesture,
  rotateAbout,
  resize,
  scaled,
  sizeRatio,
  slide,
  turn,
} from '../gestures.js';

// What a drag does to a pose: the part held stays under the cursor (moves, slides, faces, a scale by distance from a
// pivot), a turn goes at a rate.

// Looking down −Z: right is +X, up is +Y; half an Å a pixel
const frame: DragFrame = { right: Vec3.create(1, 0, 0), up: Vec3.create(0, 1, 0), perPixel: 0.5 };

const box = (rotation = Quat.identity()): GizmoPose => ({
  position: Vec3.create(0, 0, 0),
  rotation,
  scale: Vec3.create(10, 10, 10),
});

const close = (a: ArrayLike<number>, b: readonly number[]) => b.forEach((x, i) => expect(a[i]).toBeCloseTo(x, 6));

describe('gizmo gestures', () => {
  it('moves with the cursor in the screen plane (screen y is down)', () => {
    close(move(box(), frame, 20, 10).position, [10, -5, 0]);
  });

  it('slides along an axis in the screen plane at the cursor speed', () => {
    close(slide(box(), Vec3.create(1, 0, 0), frame, 20, 0).position, [10, 0, 0]);
  });

  it('slides along a foreshortened axis so the handle keeps up with the cursor', () => {
    // the axis at 60° to the screen: half its length shows, so the item goes twice as far
    const a = Vec3.normalize(Vec3(), Vec3.create(0.5, 0, Math.sqrt(3) / 2));
    expect(slide(box(), a, frame, 20, 0).position[0]).toBeCloseTo(10, 6);
  });

  it('caps a slide along an axis that points at the camera', () => {
    const p = slide(box(), Vec3.create(0, 0, 1), frame, 20, 0).position;
    expect(Vec3.magnitude(p)).toBeLessThanOrEqual(10 * MaxSlideGain);
  });

  it('slides along a local axis as the pose turns it', () => {
    // turned 90° about Z: its local +Y is world −X; dragging left slides it 10 Å that way
    const turned = box(Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2));
    close(slide(turned, Vec3.create(0, 1, 0), frame, -20, 0).position, [-10, 0, 0]);
  });

  it('drags one face, the opposite one staying', () => {
    const next = resize(box(), [{ axis: 0, sign: 1 }], frame, 8, 0);
    close(next.scale, [14, 10, 10]);
    expect(next.position[0] - next.scale[0] / 2).toBeCloseTo(-5, 6);
  });

  it('drags a face facing −Y: down the screen grows it', () => {
    const next = resize(box(), [{ axis: 1, sign: -1 }], frame, 0, 8);
    close(next.scale, [10, 14, 10]);
    expect(next.position[1] + next.scale[1] / 2).toBeCloseTo(5, 6);
  });

  it('never sizes a face below the minimum', () => {
    const next = resize(box(), [{ axis: 2, sign: 1 }], { ...frame, right: Vec3.create(0, 0, 1) }, -1000, 0);
    expect(next.scale[2]).toBeCloseTo(MinExtent, 6);
    expect(next.position[2] - next.scale[2] / 2).toBeCloseTo(-5, 6);
  });

  it('follows the rotation: a turned item sizes along its own axis', () => {
    const turned = box(Quat.setAxisAngle(Quat(), Vec3.create(0, 0, 1), Math.PI / 2));
    // its local +X is world +Y: dragging up the screen (dy < 0) grows it
    const next = resize(turned, [{ axis: 0, sign: 1 }], frame, 0, -8);
    expect(next.scale[0]).toBeCloseTo(14, 6);
    close(next.position, [0, 2, 0]);
  });

  it('turns about the screen axes, back to where it was when dragged back', () => {
    const there = turn(box(), frame, 90, 0);
    expect((Quat.getAxisAngle(Vec3(), there.rotation) * 180) / Math.PI).toBeCloseTo(45, 4);
    const back = turn(there, frame, -90, 0);
    expect(Math.abs(Quat.getAxisAngle(Vec3(), back.rotation))).toBeCloseTo(0, 4);
  });

  it('scales by distance from the pivot: outward grows on any side', () => {
    // grabbed 40 px left of the pivot: dragging further left (outward) grows it, as on the right side
    expect(sizeRatio({ ...frame, grab: [-40, 0] }, -40, 0)).toBeCloseTo(2, 6);
    expect(sizeRatio({ ...frame, grab: [40, 0] }, 40, 0)).toBeCloseTo(2, 6);
    expect(sizeRatio({ ...frame, grab: [-40, 0] }, 20, 0)).toBeCloseTo(0.5, 6);
    // grabbed at the pivot: the rate rule
    expect(sizeRatio({ ...frame, grab: [1, 0] }, 150, 0)).toBeCloseTo(2, 6);
  });

  it('scales some axes only (a box across its edges, a cylinder about its axis)', () => {
    close(scaled(box(), 1.5, [0, 2]).scale, [15, 10, 15]);
    const across = poseGesture(box(), { kind: 'scale', axes: [0, 2] }, { ...frame, grab: [30, 0] }, 15, 0)!;
    close(across.scale, [15, 10, 15]);
    close(across.position, [0, 0, 0]);
  });

  it("leaves a custom handle's gesture to the binding", () => {
    expect(poseGesture(box(), { kind: 'custom', id: 'x' }, frame, 10, 0)).toBeUndefined();
  });

  it('moves in a plane: the drag less its part along the normal', () => {
    // the plane across X: a drag right (along X) does nothing, up moves it
    close(inPlane(box(), Vec3.create(1, 0, 0), frame, 20, -10).position, [0, 5, 0]);
  });

  it('rotates about an axis facing the viewer by the angle the cursor went around it', () => {
    // Z toward the viewer; pressed 20 px right of its place, dragged a quarter turn counterclockwise (to 20 px above)
    const next = rotateAbout(box(), Vec3.create(0, 0, 1), { ...frame, grab: [20, 0] }, -20, -20);
    const axis = Vec3();
    const rad = Quat.getAxisAngle(axis, next.rotation);
    expect((rad * 180) / Math.PI).toBeCloseTo(90, 4);
    close(axis, [0, 0, 1]);
  });

  it('rotates about an axis in the screen at a rate across it', () => {
    // Y in the screen: a drag right (across it) turns it
    const next = rotateAbout(box(), Vec3.create(0, 1, 0), { ...frame, grab: [20, 0] }, 20, 0);
    expect(Math.abs((Quat.getAxisAngle(Vec3(), next.rotation) * 180) / Math.PI)).toBeCloseTo(10, 4);
  });

  const angleAbout = (pose: GizmoPose, axis: readonly number[]) => {
    const a = Vec3();
    const rad = Quat.getAxisAngle(a, pose.rotation);
    return ((rad * 180) / Math.PI) * Math.sign(Vec3.dot(a, Vec3.create(axis[0], axis[1], axis[2]))) || 0;
  };

  it('rotates by the grabbed point, facing the viewer: either side of the ring follows the cursor', () => {
    // Z toward the viewer, the ring of radius 10; a drag up 20 px (10 Å)
    const right = rotateAbout(box(), Vec3.create(0, 0, 1), { ...frame, grabPoint: Vec3.create(10, 0, 0) }, 0, -20);
    expect(angleAbout(right, [0, 0, 1])).toBeCloseTo(45, 4);
    const left = rotateAbout(box(), Vec3.create(0, 0, 1), { ...frame, grabPoint: Vec3.create(-10, 0, 0) }, 0, -20);
    expect(angleAbout(left, [0, 0, 1])).toBeCloseTo(-45, 4);
  });

  it('rotates by the grabbed point, edge-on: the near and the far side follow the cursor', () => {
    // Y in the screen, the ring in X–Z seen edge-on; a drag right 20 px (10 Å) along the radius-10 ring: 1 rad
    const near = rotateAbout(box(), Vec3.create(0, 1, 0), { ...frame, grabPoint: Vec3.create(0, 0, 10) }, 20, 0);
    expect(angleAbout(near, [0, 1, 0])).toBeCloseTo(180 / Math.PI, 4);
    const far = rotateAbout(box(), Vec3.create(0, 1, 0), { ...frame, grabPoint: Vec3.create(0, 0, -10) }, 20, 0);
    expect(angleAbout(far, [0, 1, 0])).toBeCloseTo(-180 / Math.PI, 4);
  });

  it('rotates by the grabbed point, oblique: no flip as the cursor crosses the axis', () => {
    // the axis tilted toward the viewer (cos 0.6); a drag across from a grab on the ring keeps one sign as it grows
    const axis = Vec3.normalize(Vec3(), Vec3.create(0, 0.8, 0.6));
    const grabPoint = Vec3.scale(Vec3(), Vec3.normalize(Vec3(), Vec3.cross(Vec3(), axis, Vec3.create(1, 0, 0))), 10);
    const signs = [5, 20, 40, 80].map((d) =>
      Math.sign(angleAbout(rotateAbout(box(), axis, { ...frame, grabPoint }, d, 0), axis)),
    );
    expect(new Set(signs).size).toBe(1);
  });
});
