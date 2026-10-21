/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Face, Handle } from './handle.js';

// What a drag on a handle does to an item's pose, as pure functions of the drag: pixels moved, the camera's axes and how
// many Å a pixel is at the item (`DragFrame`). Each gesture starts from the pose as it was when the drag began and the
// drag's total so far, so it never accumulates rounding. The part held stays under the cursor where a gesture can keep
// it (a move, a slide, a scale by distance from a pivot); a turn goes at a rate.

/** Where an item is, which way it is turned, how large */
export interface GizmoPose {
  readonly position: Vec3;
  readonly rotation: Quat;
  readonly scale: Vec3;
}

export interface DragFrame {
  /** The camera's right and up, unit vectors in the world */
  readonly right: Vec3;
  readonly up: Vec3;
  /** Å a CSS pixel is at the item's depth */
  readonly perPixel: number;
  /** Where the drag started from the handle's pivot on the screen (CSS px, x right, y down) */
  readonly grab?: readonly [number, number];
}

/** Degrees of turn a pixel of drag */
export const TurnPerPixel = 0.5;
/** How much faster than the cursor a slide may go along an axis that points (nearly) at the camera */
export const MaxSlideGain = 4;
/** Pixels of drag that double (or halve) a size when there is no grab to scale by */
export const DoublingPixels = 150;
/** A grab closer than this (px) to its pivot has no direction to scale along: the rate rule instead */
export const MinGrab = 6;
/** The smallest extent (Å) a drag sizes anything to */
export const MinExtent = 0.2;

/** The drag in the world: `dx` right, `dy` down on the screen */
export function worldDelta(frame: DragFrame, dx: number, dy: number): Vec3 {
  const out = Vec3.scale(Vec3(), frame.right, dx * frame.perPixel);
  return Vec3.scaleAndAdd(out, out, frame.up, -dy * frame.perPixel);
}

/** How far (Å) along the unit `axis` the drag goes, as the screen shows the axis (foreshortened, capped) */
export function alongAxis(axis: Vec3, frame: DragFrame, dx: number, dy: number): number {
  // the axis on the screen, its length 1 when it lies in the screen's plane
  const sx = Vec3.dot(axis, frame.right);
  const sy = Vec3.dot(axis, frame.up);
  const onScreen = Math.max(sx * sx + sy * sy, 1 / (MaxSlideGain * MaxSlideGain));
  return Vec3.dot(worldDelta(frame, dx, dy), axis) / onScreen;
}

/** How much farther from the grab's pivot the cursor is than the grab was; undefined grabbed (nearly) at it */
export function grabRatio(frame: DragFrame | undefined, dx: number, dy: number): number | undefined {
  const g = frame?.grab;
  const reach = g ? Math.hypot(g[0], g[1]) : 0;
  if (!g || reach < MinGrab) return undefined;
  return Math.max(Math.hypot(g[0] + dx, g[1] + dy) / reach, 0.01);
}

/** By the grab, or at a rate: twice as large every `DoublingPixels` up or right */
export function sizeRatio(frame: DragFrame | undefined, dx: number, dy: number): number {
  return grabRatio(frame, dx, dy) ?? 2 ** ((dx - dy) / DoublingPixels);
}

/** A local direction in the world, as the pose turns it */
export function worldAxis(pose: GizmoPose, local: Vec3): Vec3 {
  return Vec3.normalize(Vec3(), Vec3.transformQuat(Vec3(), local, pose.rotation));
}

/** Moved in the screen's plane */
export function move(from: GizmoPose, frame: DragFrame, dx: number, dy: number): GizmoPose {
  return { ...from, position: Vec3.add(Vec3(), from.position, worldDelta(frame, dx, dy)) };
}

/** Slid along `local` (as the pose turns it), by the drag along it as the screen shows it */
export function slide(from: GizmoPose, local: Vec3, frame: DragFrame, dx: number, dy: number): GizmoPose {
  const a = worldAxis(from, local);
  return { ...from, position: Vec3.scaleAndAdd(Vec3(), from.position, a, alongAxis(a, frame, dx, dy)) };
}

/** Moved in the plane across `localNormal` (as the pose turns it): the drag in the world, less its part along the normal */
export function inPlane(from: GizmoPose, localNormal: Vec3, frame: DragFrame, dx: number, dy: number): GizmoPose {
  const n = worldAxis(from, localNormal);
  const d = worldDelta(frame, dx, dy);
  Vec3.scaleAndAdd(d, d, n, -Vec3.dot(d, n));
  return { ...from, position: Vec3.add(Vec3(), from.position, d) };
}

/**
 * Turned about `local` (as the pose turns it) through its position. Facing the viewer, by the angle the cursor went
 * around the axis's place on the screen (`frame.grab` is the press from there): the part held stays under the cursor.
 * Lying in the screen (where that angle is unstable), by the drag across the axis, `TurnPerPixel` a pixel.
 */
export function rotateAbout(from: GizmoPose, local: Vec3, frame: DragFrame, dx: number, dy: number): GizmoPose {
  const a = worldAxis(from, local);
  const toViewer = Vec3.cross(Vec3(), frame.right, frame.up);
  const facing = Vec3.dot(a, toViewer);
  const g = frame.grab;
  let rad: number;
  if (g && Math.abs(facing) > 0.2 && Math.hypot(g[0], g[1]) >= MinGrab) {
    // on the screen with y up: the angle from the press to the cursor, counterclockwise; a turn about an axis toward
    // the viewer looks counterclockwise
    const [x0, y0] = [g[0], -g[1]];
    const [x1, y1] = [g[0] + dx, -(g[1] + dy)];
    rad = Math.atan2(x0 * y1 - y0 * x1, x0 * x1 + y0 * y1) * Math.sign(facing);
  } else {
    // the axis across the screen: a drag across it (its perpendicular on the screen, y up) turns it
    const sx = Vec3.dot(a, frame.right);
    const sy = Vec3.dot(a, frame.up);
    const len = Math.hypot(sx, sy) || 1;
    rad = (((-sy / len) * dx + (sx / len) * -dy) * TurnPerPixel * Math.PI) / 180;
  }
  const q = Quat.setAxisAngle(Quat(), a, rad);
  return { ...from, rotation: Quat.normalize(Quat(), Quat.multiply(Quat(), q, from.rotation)) };
}

/** Turned about its position: a sideways drag about the screen's up, an upward one about its right */
export function turn(from: GizmoPose, frame: DragFrame, dx: number, dy: number): GizmoPose {
  const rad = (d: number) => (d * TurnPerPixel * Math.PI) / 180;
  const about = Quat.multiply(
    Quat(),
    Quat.setAxisAngle(Quat(), frame.up, rad(dx)),
    Quat.setAxisAngle(Quat(), frame.right, rad(dy)),
  );
  return { ...from, rotation: Quat.normalize(Quat(), Quat.multiply(Quat(), about, from.rotation)) };
}

/** `axes` (all when absent) scaled by `ratio`, none below `MinExtent` */
export function scaled(from: GizmoPose, ratio: number, axes: readonly (0 | 1 | 2)[] = [0, 1, 2]): GizmoPose {
  const scale = Vec3.clone(from.scale);
  for (const a of axes) scale[a] = Math.max(MinExtent, from.scale[a] * ratio);
  return { ...from, scale };
}

/** Faces dragged along their axes, the opposite faces staying (a box's edge would move the two it joins) */
export function resize(from: GizmoPose, faces: readonly Face[], frame: DragFrame, dx: number, dy: number): GizmoPose {
  const scale = Vec3.clone(from.scale);
  const position = Vec3.clone(from.position);
  for (const { axis, sign } of faces) {
    const local = Vec3.create(0, 0, 0);
    local[axis] = sign;
    const dir = worldAxis(from, local);
    const before = scale[axis];
    scale[axis] = Math.max(MinExtent, before + alongAxis(dir, frame, dx, dy));
    Vec3.scaleAndAdd(position, position, dir, (scale[axis] - before) / 2);
  }
  return { ...from, scale, position };
}

/** The pose a handle's drag makes; undefined for a `custom` handle (the binding's own) */
export function poseGesture(
  from: GizmoPose,
  handle: Handle,
  frame: DragFrame,
  dx: number,
  dy: number,
): GizmoPose | undefined {
  switch (handle.kind) {
    case 'move':
      return move(from, frame, dx, dy);
    case 'slide':
      return slide(from, handle.axis, frame, dx, dy);
    case 'plane':
      return inPlane(from, handle.normal, frame, dx, dy);
    case 'turn':
      return turn(from, frame, dx, dy);
    case 'rotate':
      return rotateAbout(from, handle.axis, frame, dx, dy);
    case 'scale':
      return scaled(from, sizeRatio(frame, dx, dy), handle.axes);
    case 'faces':
      return resize(from, handle.faces, frame, dx, dy);
    case 'custom':
      return undefined;
  }
}
