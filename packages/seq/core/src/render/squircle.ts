/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { closePath, cubicTo, lineTo, moveTo, type PathList, Round } from './draw-lists.js';

// The squircle (_spec/README.md "Render", Shapes): one curve, defined here, that every backend draws, so Canvas2D and
// SVG draw the same corner. A corner is a quarter of the superellipse |x|⁴ + |y|⁴ = 1 as one cubic; its handles are
// where the cubic's midpoint falls on the curve (3a + 4 = 8 · 2^(−1/4)), as 0.5523 is for a circle.

/** A handle's length along its side, as a fraction of the corner's radius. */
export const SquircleHandle = (8 * 2 ** -0.25 - 4) / 3;

/**
 * The radius of a discrete mark of `h` px, by the shape rules: a third of its height (`ratio`), capped (`max`), and
 * none (square) below `min` px in either dimension.
 */
export function squircleRadius(w: number, h: number, ratio: number, max: number, min: number): number {
  if (w < min || h < min) return 0;
  return Math.min(h * ratio, max, w / 2, h / 2);
}

/** The outline of a rect with squircle corners (the ends `round` names; the others square) added to a path. */
export function squirclePath(p: PathList, x: number, y: number, w: number, h: number, r: number, round: number): void {
  const left = (round & Round.left) !== 0 ? Math.min(r, w / 2, h / 2) : 0;
  const right = (round & Round.right) !== 0 ? Math.min(r, w / 2, h / 2) : 0;
  const x1 = x + w;
  const y1 = y + h;
  const k = 1 - SquircleHandle;
  moveTo(p, x + left, y);
  lineTo(p, x1 - right, y);
  if (right > 0) cubicTo(p, x1 - right * k, y, x1, y + right * k, x1, y + right);
  lineTo(p, x1, y1 - right);
  if (right > 0) cubicTo(p, x1, y1 - right * k, x1 - right * k, y1, x1 - right, y1);
  lineTo(p, x + left, y1);
  if (left > 0) cubicTo(p, x + left * k, y1, x, y1 - left * k, x, y1 - left);
  lineTo(p, x, y + left);
  if (left > 0) cubicTo(p, x, y + left * k, x + left * k, y, x + left, y);
  closePath(p);
}
