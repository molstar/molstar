/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { describe, expect, it } from 'vitest';
import type { CutObject } from '../../cuts.js';
import { type DragFrame, MaxSlideGain, MinExtent, move, resize, slide, turn } from '../gestures.js';

// What a drag does to a cut: the handle stays under the cursor (moves, slides, faces), turns and sizes by rate.

// Looking down −Z: right is +X, up is +Y; half an Å a pixel
const frame: DragFrame = { right: Vec3.create(1, 0, 0), up: Vec3.create(0, 1, 0), perPixel: 0.5 };

const box = (): CutObject => ({
  type: 'cube',
  invert: false,
  position: Vec3.create(0, 0, 0),
  rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
  scale: Vec3.create(10, 10, 10),
  transform: undefined as never,
});

const close = (a: ArrayLike<number>, b: readonly number[]) => b.forEach((x, i) => expect(a[i]).toBeCloseTo(x, 6));

describe('gestures', () => {
  it('moves with the cursor in the screen plane (screen y is down)', () => {
    close(move(box(), frame, 20, 10).position, [10, -5, 0]);
  });

  it('slides along an axis in the screen plane at the cursor speed', () => {
    close(slide(box(), Vec3.create(1, 0, 0), frame, 20, 0).position, [10, 0, 0]);
  });

  it('slides along a foreshortened axis so the handle keeps up with the cursor', () => {
    // the axis at 60° to the screen: half its length shows, so the cut goes twice as far
    const a = Vec3.normalize(Vec3(), Vec3.create(0.5, 0, Math.sqrt(3) / 2));
    const p = slide(box(), a, frame, 20, 0).position;
    // on the screen it moved 10 Å (20 px) along X
    expect(p[0]).toBeCloseTo(10, 6);
  });

  it('caps a slide along an axis that points at the camera', () => {
    const p = slide(box(), Vec3.create(0, 0, 1), frame, 20, 0).position;
    expect(Vec3.magnitude(p)).toBeLessThanOrEqual(10 * MaxSlideGain);
  });

  it('drags one face of a box, the opposite one staying', () => {
    const next = resize(box(), 'size-x', frame, 8, 0);
    close(next.scale, [14, 10, 10]);
    close(next.position, [2, 0, 0]);
    // −X face where it was
    expect(next.position[0] - next.scale[0] / 2).toBeCloseTo(-5, 6);
  });

  it('drags the bottom (size-y faces −Y): down the screen grows it', () => {
    const next = resize(box(), 'size-y', frame, 0, 8);
    close(next.scale, [10, 14, 10]);
    expect(next.position[1] + next.scale[1] / 2).toBeCloseTo(5, 6);
  });

  it('never sizes a face below the minimum', () => {
    const next = resize(box(), 'size-z', { ...frame, right: Vec3.create(0, 0, 1) }, -1000, 0);
    expect(next.scale[2]).toBeCloseTo(MinExtent, 6);
    // and moves the center by only what it shrank
    expect(next.position[2] - next.scale[2] / 2).toBeCloseTo(-5, 6);
  });

  it('widens a cylinder about its axis: the diameter twice the drag', () => {
    const cylinder = { ...box(), type: 'cylinder' as const, scale: Vec3.create(10, 20, 10) };
    const next = resize(cylinder, 'radius', frame, 4, 0);
    close(next.scale, [14, 20, 14]);
    close(next.position, [0, 0, 0]);
  });

  it('follows the cut rotation: a turned box sizes along its own axis', () => {
    const turned = { ...box(), rotation: { axis: Vec3.create(0, 0, 1), angle: 90 } };
    // its local +X is world +Y: dragging up the screen (dy < 0) grows it
    const next = resize(turned, 'size-x', frame, 0, -8);
    expect(next.scale[0]).toBeCloseTo(14, 6);
    close(next.position, [0, 2, 0]);
  });

  it('turns about the screen axes, back to where it was when dragged back', () => {
    const there = turn(box(), frame, 90, 0);
    expect(there.rotation.angle).toBeCloseTo(45, 4);
    const back = turn(there, frame, -90, 0);
    expect(back.rotation.angle).toBeCloseTo(0, 4);
  });
});
