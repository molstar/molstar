/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { grown, PathCommand, type PathList, Round } from './draw-lists.js';

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

const Args = new Float64Array(6);

/** The outline of a rect with squircle corners (the ends `round` names; the others square) added to a path. */
export function squirclePath(p: PathList, x: number, y: number, w: number, h: number, r: number, round: number): void {
  Args[0] = x;
  Args[1] = y;
  Args[2] = w;
  Args[3] = h;
  Args[4] = r;
  Args[5] = round;
  squirclePathArgs(p, Args);
}

/**
 * `squirclePath` with its numbers in `a` (x, y, w, h, radius, round): what a draw loop calls, since a number passed to
 * a call that is not inlined is boxed (the polar warp traces every rect through it: nothing allocated once warm).
 */
export function squirclePathArgs(p: PathList, a: Float64Array): void {
  const x = a[0];
  const y = a[1];
  const w = a[2];
  const h = a[3];
  const r = Math.min(a[4], w / 2, h / 2);
  const round = a[5];
  // each corner round where its side's end is (`Round.left`, `right`) and its edge is not cut (`cutTop`, `cutBottom`)
  const left = (round & Round.left) !== 0;
  const right = (round & Round.right) !== 0;
  const top = (round & Round.cutTop) === 0;
  const bottom = (round & Round.cutBottom) === 0;
  const tl = left && top ? r : 0;
  const tr = right && top ? r : 0;
  const br = right && bottom ? r : 0;
  const bl = left && bottom ? r : 0;
  const x1 = x + w;
  const y1 = y + h;
  const k = 1 - SquircleHandle;
  // the outline's commands and points written straight into the path's arrays: a point passed to a call that is not
  // inlined is boxed (the polar warp traces every rect through here)
  const P = Pts;
  let n = 0;
  let c = 0;
  Cmds[c++] = PathCommand.move;
  P[n++] = x + tl;
  P[n++] = y;
  Cmds[c++] = PathCommand.line;
  P[n++] = x1 - tr;
  P[n++] = y;
  if (tr > 0) {
    Cmds[c++] = PathCommand.cubic;
    P[n++] = x1 - tr * k;
    P[n++] = y;
    P[n++] = x1;
    P[n++] = y + tr * k;
    P[n++] = x1;
    P[n++] = y + tr;
  }
  Cmds[c++] = PathCommand.line;
  P[n++] = x1;
  P[n++] = y1 - br;
  if (br > 0) {
    Cmds[c++] = PathCommand.cubic;
    P[n++] = x1;
    P[n++] = y1 - br * k;
    P[n++] = x1 - br * k;
    P[n++] = y1;
    P[n++] = x1 - br;
    P[n++] = y1;
  }
  Cmds[c++] = PathCommand.line;
  P[n++] = x + bl;
  P[n++] = y1;
  if (bl > 0) {
    Cmds[c++] = PathCommand.cubic;
    P[n++] = x + bl * k;
    P[n++] = y1;
    P[n++] = x;
    P[n++] = y1 - bl * k;
    P[n++] = x;
    P[n++] = y1 - bl;
  }
  Cmds[c++] = PathCommand.line;
  P[n++] = x;
  P[n++] = y + tl;
  if (tl > 0) {
    Cmds[c++] = PathCommand.cubic;
    P[n++] = x;
    P[n++] = y + tl * k;
    P[n++] = x + tl * k;
    P[n++] = y;
    P[n++] = x + tl;
    P[n++] = y;
  }
  Cmds[c++] = PathCommand.close;
  append(p, c, n);
}

/** A squircle's outline at most: 10 commands (a move, four lines, four corners, a close), 17 points. */
const Cmds = new Uint8Array(10);
const Pts = new Float64Array(34);

/** Appends the outline's `c` commands and `n` numbers (x, y pairs) to a path. */
function append(p: PathList, c: number, n: number): void {
  p.commands = grown(p.commands, p.commandCount + c);
  for (let i = 0; i < c; i++) p.commands[p.commandCount++] = Cmds[i];
  p.points = grown(p.points, 2 * p.pointCount + n);
  const at = 2 * p.pointCount;
  for (let i = 0; i < n; i++) p.points[at + i] = Pts[i];
  p.pointCount += n >> 1;
}
