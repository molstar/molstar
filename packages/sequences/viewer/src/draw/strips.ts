/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Strip lanes (seqspec:viewer#the-text-view), part of the renderer API: a `values` or `color` layer alone on its track
// in the text view is a strip under each line, `wrappedStripHeight` high after `wrappedTrackGap` (the band
// `stripLineHeights` gives), a color per position from the layer's theme; beside the letters (a taller band) the same
// colors are the letters' cells. Drawn by run of one color within a group, cut at the gaps between groups; the hovered
// position's strip cell outlined (seqspec:interaction#hover). Allocation-free: module scratch.

import { NO_COLOR, type RectBatch } from './backend.js';
import { cellBand } from './cells.js';
import { outlineItem } from './hover.js';
import type { DrawView } from './view.js';

/** A strip's colors: per block index (its shown index beside, for a source read by shown position: a stream's). */
export interface StripColors {
  /** 0xRRGGBB, a skin token, or `NO_COLOR` (none: nothing drawn there). */
  color(block: number, index: number, shown: number): number;
  /** The alpha byte of the color `color` gave last. */
  alpha(block: number, index: number): number;
}

/**
 * Wrapped: whether the view's band is a strip lane (no taller than `wrappedTrackGap` + `wrappedStripHeight`), with
 * `cellBand` its rows (device px) at the band's bottom; else false, `cellBand` the band's cells (2 px off both sides of a
 * band over 8 px), behind the letters the band holds.
 */
export function stripBand(v: DrawView): boolean {
  const strip = v.height <= v.trackGap + v.stripHeight + 0.001;
  const pad = strip || v.height <= 8 ? 0 : 2;
  const y1 = Math.round((v.top + v.height - pad) * v.dpr) | 0;
  const y0 = strip
    ? Math.min(y1 - 1, Math.round((v.top + v.height - v.stripHeight) * v.dpr) | 0)
    : Math.round((v.top + pad) * v.dpr) | 0;
  cellBand[0] = y0;
  cellBand[1] = y1 - y0;
  cellBand[2] = 0;
  return strip;
}

/**
 * Wrapped: block b's shown [lo, hi) on the view's line in `cellBand`'s rows, into `rects`: one rect per run of one color
 * and alpha, a run cut where the next position does not touch it (within a device px: a gap between groups); a position with no color
 * draws nothing.
 */
export function drawStrip(v: DrawView, b: number, lo: number, hi: number, src: StripColors, rects: RectBatch): void {
  const ax = v.axis;
  if (!ax || lo >= hi) return;
  const y = cellBand[0];
  const h = cellBand[1];
  let x0 = 0;
  let x1 = 0;
  let c = NO_COLOR;
  let a = 0;
  for (let k = lo; k < hi; k++) {
    const i = ax.indexOf(b, k as never);
    const color = src.color(b, i, k);
    const xa = v.xd(b, k);
    const xb = v.xe(b, k);
    if (color === NO_COLOR) {
      if (c !== NO_COLOR) rects.push(x0, y, Math.max(1, x1 - x0), h);
      c = NO_COLOR;
      continue;
    }
    const alpha = src.alpha(b, i);
    if (c !== NO_COLOR && color === c && alpha === a && xa <= x1 + 1) {
      x1 = xb;
      continue;
    }
    if (c !== NO_COLOR) rects.push(x0, y, Math.max(1, x1 - x0), h);
    c = color;
    a = alpha;
    x0 = xa;
    x1 = xb;
    rects.begin(c, a);
  }
  if (c !== NO_COLOR) rects.push(x0, y, Math.max(1, x1 - x0), h);
}

const el = new Int32Array(4);
/**
 * Wrapped: the hovered position (the layer's own hit, `marks.element`) on a strip lane of this line, its strip cell
 * outlined inside in its color's hover ink; nothing beside the letters (their own hover shows it) or where it has no
 * color.
 */
export function outlineStrip(v: DrawView, src: StripColors): void {
  const ax = v.axis;
  if (!ax || !v.marks.element(v.layer, el) || el[0] !== v.onlyBlock || !stripBand(v)) return;
  const b = el[0];
  const k = v.shownOf(b, el[1]);
  if (k < v.lineStart || k >= v.lineStop) return;
  const color = src.color(b, el[1], k);
  if (color === NO_COLOR) return;
  const x = v.xd(b, k);
  outlineItem(v, color, x, cellBand[0], Math.max(1, v.xe(b, k) - x), cellBand[1]);
}
