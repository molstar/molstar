/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The placement helpers (seqspec:viewer#extension-points): a draw that places through them draws in the linear, wrapped
// and circular geometries alike (rings draw the axis unrolled, mapped by the composer's `RingBackend`). Runs are in the
// block coordinate system and go through the axis's run converters (seqspec:model#pieces-and-layout-maps), so a
// reversed piece, a filter or a cyclic origin costs a draw nothing. Outputs are device px (x edges rounded, as
// rectangles take them); `y` inputs are the view's CSS px (`v.top` and down).
// Allocation-free: module scratch, grown only for a run that crosses more filter runs than any before.

import { AxisGeometry, type AxisView, type ShownIndex } from '../axis.js';
import type { RectBatch } from './backend.js';
import type { DrawView } from './view.js';

let shown = new Int32Array(32);
const piece = new Float64Array(64);

/**
 * Block b's run [start, stop) as drawn here, its pieces into `out` (at most `out.length / 2`), their count: linear, the
 * x edges of each shown run in view; wrapped, those on the view's line.
 */
export function placeRun(v: DrawView, b: number, start: number, stop: number, out: Float64Array): number {
  if ((v.onlyBlock >= 0 && b !== v.onlyBlock) || stop <= start) return 0;
  let n = v.runToShown(b, start, stop, shown);
  if (2 * n > shown.length) {
    shown = new Int32Array(4 * n);
    n = v.runToShown(b, start, stop, shown);
  }
  const g = v.geometry;
  const lo = g === AxisGeometry.wrapped ? v.lineStart : v.visLo[b];
  const hi = g === AxisGeometry.wrapped ? v.lineStop : v.visHi[b];
  let m = 0;
  for (let r = 0; r < n && 2 * m + 1 < out.length; r++) {
    const k0 = Math.max(lo, shown[2 * r]);
    const k1 = Math.min(hi, shown[2 * r + 1]);
    if (k1 <= k0) continue;
    // The device edges as int32 (`xd`, `xe`: no double crosses a call, which the guards caught allocating a boxed one
    // per edge on a path called a few times a frame, as the overlay's skeletons, and wrapped per residue placed).
    out[2 * m] = v.xd(b, k0);
    // Wrapped, the last cell's right edge (the next position may start a group past a gap).
    out[2 * m + 1] = g === AxisGeometry.wrapped ? v.xe(b, k1 - 1) : v.xd(b, k1);
    m++;
  }
  return m;
}

/** A wrapped position's advance: its prepared run's (the font's cell × its span). */
const advance = (v: DrawView, b: number) => (v.runs[b] ? v.runs[b].advance : v.cell * v.span[b]);

/** Block b's run [start, stop) between y0 and y1 (CSS px): a rect per piece (on rings, mapped to a sector). */
export function band(v: DrawView, b: number, start: number, stop: number, y0: number, y1: number, rects: RectBatch) {
  const n = placeRun(v, b, start, stop, piece);
  const top = Math.round(y0 * v.dpr);
  const h = Math.max(1, Math.round(y1 * v.dpr) - top);
  for (let p = 0; p < n; p++) rects.push(piece[2 * p], top, Math.max(1, piece[2 * p + 1] - piece[2 * p]), h);
}

/**
 * The point to draw a glyph or text at for block b's index i at the view's y (CSS px), device px into `out[0..1]`:
 * the center of its cell. False, `out` untouched: not shown here.
 */
export function mark(v: DrawView, b: number, i: number, y: number, out: Float64Array): boolean {
  if (v.onlyBlock >= 0 && b !== v.onlyBlock) return false;
  const k = v.shownOf(b, i);
  const g = v.geometry;
  if (k < 0) return false;
  if (g === AxisGeometry.wrapped ? k < v.lineStart || k >= v.lineStop : k < v.visLo[b] || k >= v.visHi[b]) return false;
  out[0] =
    (v.x(b, k) + (g === AxisGeometry.wrapped ? advance(v, b) : v.x(b, (k + 1) as ShownIndex) - v.x(b, k)) / 2) * v.dpr;
  out[1] = y * v.dpr;
  return true;
}

/**
 * The shown edge (a shown index: that cell's left edge) of the boundary before block b's index `i` (seqspec:format#placement):
 * its position's leading edge in reading order (a reversed piece's right). A boundary after the system's last position
 * (a C-terminal insertion), or before a position the block's filter hides, sits on the trailing edge of the position
 * before it; -1: not in this block (a junction between two pieces of one unit is drawn once, by the piece holding the
 * position after it).
 */
export function boundaryEdge(v: AxisView, b: number, i: number): number {
  const q = v.shownOf(b, i);
  if (q >= 0) return v.dir[b] < 0 ? q + 1 : q;
  const ax = v.axis;
  if (!ax || i <= 0 || (i !== ax.base[b] + ax.length[b] && ax.positions[b] === null)) return -1;
  const p = v.shownOf(b, i - 1);
  return p < 0 ? -1 : v.dir[b] < 0 ? p : p + 1;
}
