/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { type CutObject, coneAngle, coneScale, rotated } from '../cuts.js';

// What a drag on a gizmo's handle does to its cut (clipping.md §4), as pure functions of the drag: pixels moved,
// the camera's axes and how many Å a pixel is at the cut. Each gesture starts from the cut as it was when the drag
// began (`from`) and the drag's total so far, so it never accumulates rounding.

/** A face of a box (or a cylinder's end): its local axis and which way it faces */
export interface Face {
  readonly axis: 0 | 1 | 2;
  readonly sign: 1 | -1;
}

/**
 * A handle of a gizmo: what dragging it does. `size` scales the whole cut (a sphere's ring, a box's corner); `faces`
 * drags faces along their axes, the opposite ones staying (a box's edge: the two it joins; a cylinder's end); `radius`
 * widens a cylinder and `angle` a cone, about `pivot` (in the cut's frame: a rim's center).
 */
export type Handle =
  | { readonly kind: 'move' }
  | { readonly kind: 'slide' }
  | { readonly kind: 'turn' }
  | { readonly kind: 'size' }
  | { readonly kind: 'faces'; readonly faces: readonly Face[] }
  | { readonly kind: 'radius'; readonly pivot: Vec3 }
  | { readonly kind: 'angle'; readonly pivot: Vec3 };

/** The smallest a cut's extent gets by dragging (Å) */
export const MinExtent = 0.2;

export interface DragFrame {
  /** The camera's right and up, unit vectors in the world */
  readonly right: Vec3;
  readonly up: Vec3;
  /** Å a pixel is at the cut's depth */
  readonly perPixel: number;
  /**
   * Where the drag started from the cut's center on the screen (CSS px, x right, y down): a size drag scales by how
   * much farther from the center the cursor is, so the part grabbed stays under it on any side of the cut
   */
  readonly grab?: readonly [number, number];
}

/** A grab closer than this (px) to the center has no direction to scale along: the rate rule instead */
export const MinGrab = 6;

/** Degrees of turn a pixel of drag */
export const TurnPerPixel = 0.5;
/** How much faster than the cursor a slide may go along an axis that points (nearly) at the camera */
export const MaxSlideGain = 4;
/** Pixels of drag that double (or halve) a size */
export const DoublingPixels = 150;

/** The drag in the world: `dx` right, `dy` down on the screen */
export function worldDelta(frame: DragFrame, dx: number, dy: number): Vec3 {
  const out = Vec3.scale(Vec3(), frame.right, dx * frame.perPixel);
  return Vec3.scaleAndAdd(out, out, frame.up, -dy * frame.perPixel);
}

/** The cut moved in the screen's plane */
export function move(from: CutObject, frame: DragFrame, dx: number, dy: number): CutObject {
  return { ...from, position: Vec3.add(Vec3(), from.position, worldDelta(frame, dx, dy)) };
}

/**
 * The cut slid along `axis` (a plane along its normal) so that the handle stays under the cursor: by the drag along the
 * axis as the screen shows it, foreshortened (an axis pointing at the camera barely moves: capped at `MaxSlideGain`)
 */
export function slide(from: CutObject, axis: Vec3, frame: DragFrame, dx: number, dy: number): CutObject {
  const a = Vec3.normalize(Vec3(), axis);
  return { ...from, position: Vec3.scaleAndAdd(Vec3(), from.position, a, alongAxis(a, frame, dx, dy)) };
}

/** How far (Å) along the unit `axis` the drag goes, as the screen shows the axis (foreshortened, capped) */
export function alongAxis(a: Vec3, frame: DragFrame, dx: number, dy: number): number {
  // the axis on the screen, its length 1 when it lies in the screen's plane
  const sx = Vec3.dot(a, frame.right);
  const sy = Vec3.dot(a, frame.up);
  const onScreen = Math.max(sx * sx + sy * sy, 1 / (MaxSlideGain * MaxSlideGain));
  return Vec3.dot(worldDelta(frame, dx, dy), a) / onScreen;
}

/** Faces dragged along their axes, the opposite faces staying: a box's edge moves the two it joins */
export function resize(from: CutObject, faces: readonly Face[], frame: DragFrame, dx: number, dy: number): CutObject {
  const scale = Vec3.clone(from.scale);
  const position = Vec3.clone(from.position);
  for (const { axis, sign } of faces) {
    const local = Vec3.create(0, 0, 0);
    local[axis] = sign;
    const dir = rotated(from.rotation, local);
    const before = scale[axis];
    scale[axis] = Math.max(MinExtent, before + alongAxis(dir, frame, dx, dy));
    Vec3.scaleAndAdd(position, position, dir, (scale[axis] - before) / 2);
  }
  return { ...from, scale, position };
}

/** How much farther from the grab's center the cursor is than the grab was; undefined grabbed (nearly) at it */
export function grabRatio(frame: DragFrame | undefined, dx: number, dy: number): number | undefined {
  const g = frame?.grab;
  const reach = g ? Math.hypot(g[0], g[1]) : 0;
  if (!g || reach < MinGrab) return undefined;
  return Math.max(Math.hypot(g[0] + dx, g[1] + dy) / reach, 0.01);
}

/** A cylinder widened about its axis: its rim stays under the cursor */
export function radius(from: CutObject, frame: DragFrame, dx: number, dy: number): CutObject {
  const ratio = grabRatio(frame, dx, dy) ?? 2 ** ((dx - dy) / DoublingPixels);
  const d = Math.max(MinExtent, from.scale[0] * ratio);
  return { ...from, scale: Vec3.create(d, from.scale[1], d) };
}

/** A cone widened or narrowed: its base rim stays under the cursor (the rim's radius goes as the angle's tangent) */
export function angle(from: CutObject, frame: DragFrame, dx: number, dy: number): CutObject {
  const ratio = grabRatio(frame, dx, dy) ?? 2 ** ((dx - dy) / DoublingPixels);
  const tan = Math.tan((coneAngle(from.scale) * Math.PI) / 180) * ratio;
  const a = Math.max(1, Math.min(89, (Math.atan(tan) * 180) / Math.PI));
  return { ...from, scale: coneScale(a) };
}

/** The cut turned about its center: a sideways drag about the screen's up, an upward one about its right */
export function turn(from: CutObject, frame: DragFrame, dx: number, dy: number): CutObject {
  const rad = (d: number) => (d * TurnPerPixel * Math.PI) / 180;
  // a stored axis need not be unit (set from code): a quaternion needs one
  const axis = Vec3.normalize(Vec3(), from.rotation.axis);
  const q = Quat.setAxisAngle(Quat(), axis, (from.rotation.angle * Math.PI) / 180);
  const about = Quat.multiply(
    Quat(),
    Quat.setAxisAngle(Quat(), frame.up, rad(dx)),
    Quat.setAxisAngle(Quat(), frame.right, rad(dy)),
  );
  const next = Quat.multiply(Quat(), about, q);
  return { ...from, rotation: rotationOf(next) };
}

/**
 * The cut sized whole. Grabbed off its center (`frame.grab`): by the cursor's distance from the center over the
 * grab's, so the ring or corner held stays under the cursor whichever side it is on; else at a rate, twice as large
 * every `DoublingPixels` up or right.
 */
export function size(from: CutObject, dx: number, dy: number, frame?: DragFrame): CutObject {
  if (from.type === 'plane') return from;
  if (from.type === 'infiniteCone') return frame ? angle(from, frame, dx, dy) : from;
  const ratio = grabRatio(frame, dx, dy) ?? 2 ** ((dx - dy) / DoublingPixels);
  return { ...from, scale: Vec3.scale(Vec3(), from.scale, ratio) };
}

/** A unit quaternion as Mol*'s clip rotation (axis, degrees) */
export function rotationOf(q: Quat): Clip.Rotation {
  const axis = Vec3();
  const rad = Quat.getAxisAngle(axis, q);
  if (!Number.isFinite(rad) || Vec3.magnitude(axis) < 1e-9) return { axis: Vec3.create(1, 0, 0), angle: 0 };
  return { axis: Vec3.normalize(axis, axis), angle: (rad * 180) / Math.PI };
}

/** The gesture of a handle */
export function gesture(handle: Handle, from: CutObject, axis: Vec3, frame: DragFrame, dx: number, dy: number) {
  switch (handle.kind) {
    case 'move':
      return move(from, frame, dx, dy);
    case 'slide':
      return slide(from, axis, frame, dx, dy);
    case 'turn':
      return turn(from, frame, dx, dy);
    case 'size':
      return size(from, dx, dy, frame);
    case 'faces':
      return resize(from, handle.faces, frame, dx, dy);
    case 'radius':
      return radius(from, frame, dx, dy);
    case 'angle':
      return angle(from, frame, dx, dy);
  }
}
