/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NoColor, type PassWriter } from './draw-lists.js';

// The square's window over one section of blocks with a switch (_spec/README.md "Render", Navigators, The square, The
// window as the blocks shown): what the section shows of the lattice, a run on the diagonal (one box, its gestures the
// viewer's) or scattered (each shown block's diagonal cell and the cross cells between them, outlined cell by cell, no
// grips). Pure: the viewer's gestures and the export read the same shape.

/** A run of lattice blocks `[i0, i1)` on the diagonal, or the blocks are scattered (`i0` -1). */
export interface SquareRun {
  i0: number;
  i1: number;
}

/**
 * The shown blocks' shape: `out` set to the run `[i0, i1)` when they are one contiguous run (every block: the whole
 * lattice), else `i0 = i1 = -1` (scattered, or none shown). Returns whether they form a run.
 */
export function shownRun(shown: readonly boolean[], out: SquareRun): boolean {
  let i0 = -1;
  let i1 = -1;
  for (let k = 0; k < shown.length; k++) {
    if (!shown[k]) continue;
    if (i0 < 0) i0 = k;
    else if (i1 !== k) {
      out.i0 = out.i1 = -1;
      return false;
    }
    i1 = k + 1;
  }
  out.i0 = i0;
  out.i1 = i0 < 0 ? -1 : i1;
  return i0 >= 0;
}

/** The lattice's edges as px along one axis: `edge(k)` the px of block `k`'s start (`k = n` the end). */
export type SquareEdge = (k: number, horizontal: boolean) => number;

/**
 * The window over the blocks shown: a run as one box (its corners' grips when `style[2]` > 0), else each cell of two
 * shown blocks outlined (k² cells, within `maxCells`); in `out`'s pass, with the ink `style[0]` and the stroke
 * `style[1]` (an array: its numbers are never boxed).
 */
export function drawShownSet(
  out: PassWriter,
  shown: readonly boolean[],
  run: SquareRun,
  edge: SquareEdge,
  style: Float64Array,
  maxCells: number,
): void {
  const ink = style[0];
  const w = style[1];
  const grips = style[2];
  const g = out.args;
  if (run.i0 >= 0) {
    const x0 = edge(run.i0, true);
    const x1 = edge(run.i1, true);
    const y0 = edge(run.i0, false);
    const y1 = edge(run.i1, false);
    rectStroke(g, x0, y0, x1 - x0, y1 - y0, ink, w);
    out.rectArgs();
    // the grips: a filled square on the two diagonal corners (the run's start and end)
    if (grips > 0) {
      for (let c = 0; c < 2; c++) {
        const cx = c === 0 ? x0 : x1;
        const cy = c === 0 ? y0 : y1;
        g[0] = cx - grips / 2;
        g[1] = cy - grips / 2;
        g[2] = grips;
        g[3] = grips;
        g[4] = ink;
        g[5] = NoColor;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
    }
    return;
  }
  let n = 0;
  for (let k = 0; k < shown.length; k++) if (shown[k]) n++;
  if (n * n > maxCells) return;
  for (let i = 0; i < shown.length; i++) {
    if (!shown[i]) continue;
    const x0 = edge(i, true);
    const x1 = edge(i + 1, true);
    for (let j = 0; j < shown.length; j++) {
      if (!shown[j]) continue;
      const y0 = edge(j, false);
      const y1 = edge(j + 1, false);
      rectStroke(g, x0, y0, x1 - x0, y1 - y0, ink, w);
      out.rectArgs();
    }
  }
}

function rectStroke(g: Float64Array, x: number, y: number, w: number, h: number, ink: number, sw: number) {
  g[0] = x;
  g[1] = y;
  g[2] = Math.max(sw, w);
  g[3] = Math.max(sw, h);
  g[4] = NoColor;
  g[5] = ink;
  g[6] = sw;
  g[7] = 0;
  g[8] = 0;
}

/** A grip of the run's box under a point: its body (move), its start (left or top edge) or end (right or bottom). */
export type SquareGrip = 'body' | 'start' | 'end';

/** Which edge a grip is on: along x (a left or right edge), along y (top or bottom), or a corner (both). */
export type SquareGripAxis = 'x' | 'y' | 'both';

/** A grip under a point: what it moves and along which axis (the cursor says which). */
export interface SquareGripHit {
  grip: SquareGrip;
  axis: SquareGripAxis;
}

/**
 * The grip of the box `[x0, x1] × [y0, y1]` under a point, within `slop` px of an edge (at most a quarter of the box,
 * so a small box keeps a body); undefined outside it. The off-diagonal corners (top right, bottom left) go to the
 * nearer edge.
 */
export function gripOf(
  x: number,
  y: number,
  x0: number,
  x1: number,
  y0: number,
  y1: number,
  slop: number,
  out: SquareGripHit,
): boolean {
  if (x < x0 - slop || x > x1 + slop || y < y0 - slop || y > y1 + slop) return false;
  const sx = Math.min(slop, (x1 - x0) / 4);
  const sy = Math.min(slop, (y1 - y0) / 4);
  const left = x <= x0 + sx;
  const right = x >= x1 - sx;
  const top = y <= y0 + sy;
  const bottom = y >= y1 - sy;
  if ((left && top) || (right && bottom)) {
    out.grip = left ? 'start' : 'end';
    out.axis = 'both';
  } else if ((right && top) || (left && bottom)) {
    // an off-diagonal corner: the nearer edge's
    const dx = right ? x1 - x : x - x0;
    const dy = top ? y - y0 : y1 - y;
    if (dx <= dy) {
      out.grip = right ? 'end' : 'start';
      out.axis = 'x';
    } else {
      out.grip = top ? 'start' : 'end';
      out.axis = 'y';
    }
  } else if (left || right) {
    out.grip = left ? 'start' : 'end';
    out.axis = 'x';
  } else if (top || bottom) {
    out.grip = top ? 'start' : 'end';
    out.axis = 'y';
  } else {
    if (x < x0 || x > x1 || y < y0 || y > y1) return false;
    out.grip = 'body';
    out.axis = 'both';
  }
  return true;
}
