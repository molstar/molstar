/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Quat, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { type CutObject, coneAngle, coneScale, rotated } from '../cuts.js';

// What a drag on a gizmo's handle does to its cut (_spec/README.md §4), as pure functions of the drag: pixels moved,
// the camera's axes and how many Å a pixel is at the cut. Each gesture starts from the cut as it was when the drag
// began (`from`) and the drag's total so far, so it never accumulates rounding.

/**
 * A handle of a gizmo: what dragging it does. `size` scales the whole cut; `size-x`, `size-y`, `size-z` drag one face
 * of a box (or a cylinder's bottom) along its axis, the opposite face staying; `radius` widens a cylinder about its axis
 */
export type Handle = 'move' | 'slide' | 'turn' | 'size' | 'size-x' | 'size-y' | 'size-z' | 'radius';

/** The face a sizing handle drags: the local axis and which way it faces (`size-y` is the bottom: +Y has the arrow) */
export const SizeFaces = {
  'size-x': { axis: 0, sign: 1 },
  'size-y': { axis: 1, sign: -1 },
  'size-z': { axis: 2, sign: 1 },
  radius: { axis: 0, sign: 1 },
} as const;

/** The smallest a cut's extent gets by dragging (Å) */
export const MinExtent = 0.2;

export interface DragFrame {
  /** The camera's right and up, unit vectors in the world */
  readonly right: Vec3;
  readonly up: Vec3;
  /** Å a pixel is at the cut's depth */
  readonly perPixel: number;
}

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

/** A face of a box (or a cylinder's bottom) dragged along its axis, the opposite face staying; a cylinder widened */
export function resize(from: CutObject, handle: keyof typeof SizeFaces, frame: DragFrame, dx: number, dy: number) {
  const { axis, sign } = SizeFaces[handle];
  const local = Vec3.create(0, 0, 0);
  local[axis] = sign;
  const dir = rotated(from.rotation, local);
  const along = alongAxis(dir, frame, dx, dy);
  const scale = Vec3.clone(from.scale);
  if (handle === 'radius') {
    // about the axis: the diameter grows twice the drag
    const d = Math.max(MinExtent, from.scale[0] + 2 * along);
    scale[0] = d;
    scale[2] = d;
    return { ...from, scale };
  }
  scale[axis] = Math.max(MinExtent, from.scale[axis] + along);
  const grown = scale[axis] - from.scale[axis];
  return { ...from, scale, position: Vec3.scaleAndAdd(Vec3(), from.position, dir, grown / 2) };
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

/** The cut sized: twice as large every `DoublingPixels` up or right; a cone's angle widened instead */
export function size(from: CutObject, dx: number, dy: number): CutObject {
  const steps = (dx - dy) / DoublingPixels;
  if (from.type === 'infiniteCone') {
    const angle = Math.max(1, Math.min(89, coneAngle(from.scale) + steps * 15));
    return { ...from, scale: coneScale(angle) };
  }
  if (from.type === 'plane') return from;
  return { ...from, scale: Vec3.scale(Vec3(), from.scale, 2 ** steps) };
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
  switch (handle) {
    case 'move':
      return move(from, frame, dx, dy);
    case 'slide':
      return slide(from, axis, frame, dx, dy);
    case 'turn':
      return turn(from, frame, dx, dy);
    case 'size':
      return size(from, dx, dy);
    default:
      return resize(from, handle, frame, dx, dy);
  }
}
