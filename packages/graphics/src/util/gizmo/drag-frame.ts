/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Camera } from '../../canvas3d/camera.js';
import type { Canvas3D } from '../../canvas3d/canvas3d.js';
import { Viewport } from '../../canvas3d/camera/util.js';
import type { DragFrame } from './gestures.js';

// A drag's frame from the canvas at the press: the camera's axes, Å per CSS pixel at the item's depth, where the press
// is from a pivot. Positions are CSS pixels on the canvas (`InputObserver`'s `x`, `y`; its `dx`, `dy` are fractions of
// the canvas, and its first drag event comes a few pixels after the press).

const tmpFrame = Viewport();

/**
 * The size in CSS pixels of the frame the camera's field of view spans (`Camera.getFrame`): the viewport's, or the
 * whole image's where the view is cropped (a tile of a larger image)
 */
function frameCssSize(canvas: Canvas3D): [number, number] {
  const ratio = canvas.input.pixelRatio || 1;
  const { width, height } = Camera.getFrame(tmpFrame, canvas.camera);
  return [width / ratio, height / ratio];
}

/** The canvas element's size in CSS pixels, what the input's moves are fractions of */
function elementCssSize(canvas: Canvas3D): [number, number] {
  const { input } = canvas;
  const ratio = input.pixelRatio || 1;
  return [input.width / ratio, input.height / ratio];
}

/** Where a drag was pressed: its first event's position less that event's move (fractions of the canvas element) */
export function pressPoint(canvas: Canvas3D, x: number, y: number, fx: number, fy: number): [number, number] {
  const [w, h] = elementCssSize(canvas);
  return [x - fx * w, y - fy * h];
}

/** The camera's right and up, unit vectors in the world */
export function screenAxes(canvas: Canvas3D): { right: Vec3; up: Vec3 } {
  const s = canvas.camera.state;
  const dir = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), s.target, s.position));
  const right = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), dir, s.up));
  const up = Vec3.normalize(Vec3(), Vec3.cross(Vec3(), right, dir));
  return { right, up };
}

/** Å a CSS pixel is at `at`: the height the camera sees there over the height of the frame it spans */
export function perPixel(canvas: Canvas3D, at: Vec3): number {
  const s = canvas.camera.state;
  const toTarget = Vec3.normalize(Vec3(), Vec3.sub(Vec3(), s.target, s.position));
  const depth =
    s.mode === 'orthographic'
      ? Vec3.distance(s.position, s.target)
      : Math.max(1e-3, Vec3.dot(Vec3.sub(Vec3(), at, s.position), toTarget));
  return (2 * Math.tan(s.fov / 2) * depth) / frameCssSize(canvas)[1];
}

/** The press (CSS px) from `pivot` as the screen shows it: x right, y down */
export function grabFrom(canvas: Canvas3D, pivot: Vec3, x: number, y: number): [number, number] {
  const ratio = canvas.input.pixelRatio || 1;
  const p = canvas.camera.project(Vec4(), pivot);
  return [x - p[0] / ratio, y - p[1] / ratio];
}

/**
 * Radians the trackball turns the camera a CSS pixel: its `getRotateFactor` (the speed, the pixel ratio, the element's
 * aspect) over its mouse-on-circle coordinates, a pixel being 2 / the width of the frame it measures by
 * (`controls.viewport`, `Camera.getFrame`); the element's height where the frame is the whole element
 */
export function trackballTurnRate(canvas: Canvas3D): number {
  const [ew, eh] = elementCssSize(canvas);
  const aspect = ew / eh || 1;
  return (2 * canvas.props.trackball.rotateSpeed * aspect) / frameCssSize(canvas)[0];
}

/** A drag's frame: the camera's axes, Å a pixel at `at`, the press from `pivot`, the trackball's rate */
export function dragFrame(canvas: Canvas3D, at: Vec3, pivot: Vec3, x0: number, y0: number): DragFrame {
  return {
    ...screenAxes(canvas),
    perPixel: perPixel(canvas, at),
    grab: grabFrom(canvas, pivot, x0, y0),
    turnRate: trackballTurnRate(canvas),
  };
}
