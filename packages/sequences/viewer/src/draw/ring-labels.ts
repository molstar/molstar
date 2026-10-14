/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Item labels on rings (seqspec:viewer#the-render-backend, CG2): a circular section's labeled items (a renderer's
// `ringLabels`: `ranges` with `labels`) are set outside the outermost ring (past an outermost ruler's numbers), each
// with a leader: from its item's ring out along its own angle to the rings' outer edge, across the numbers to the
// label's own angle, then to the label, so a cluster's leaders leave the ring apart and fan out over the numbers, never
// through one point. Placed headlessly when the layout or the view changes: sorted by angle and spread to a line apart
// per side (right and left of the center), one sweep each way, O(n log n); labels a side cannot hold are dropped by
// priority (the layer's order, then the longer item). Drawn per ring frame from the kept arrays, upright, on the real
// backend (not through the ring map): nothing is allocated per frame.

import { fitIn, fitted, type TextManagerApi } from '../text.js';
import type { RenderBackend } from './backend.js';
import { PathBatch, PathMode, PathVerb, SkinToken, TextAlign, TextBaseline, TextBatch } from './backend.js';

/**
 * What a renderer lists for the ring's labels: an item's extent on block `block`, its text, and its annotation and row
 * (an item's runs across a cyclic origin are labeled once, at the longer).
 */
export interface RingLabelSink {
  add(block: number, start: number, stop: number, text: string, annotation: number, row: number): void;
}

/**
 * A candidate: its angle (radians), its ring's outer radius (CSS px), its text and priority (lower first); its item
 * (layer, annotation, row) and a point of its track on the unrolled axis (`ux`, `uy`, CSS px), which a hit on the
 * label is hit as (CG3).
 */
export interface RingLabelInput {
  readonly angle: number;
  readonly r0: number;
  readonly text: string;
  readonly priority: number;
  readonly layer?: number;
  readonly annotation?: number;
  readonly row?: number;
  readonly ux?: number;
  readonly uy?: number;
}

/**
 * The ring's frame for placement, CSS px: the center, the rings' outer edge (`edge`, where leaders leave them), the
 * labels' radius, the box they stay in, a line's height.
 */
export interface RingLabelFrame {
  readonly cx: number;
  readonly cy: number;
  readonly edge: number;
  readonly r: number;
  readonly left: number;
  readonly right: number;
  readonly top: number;
  readonly bottom: number;
  readonly line: number;
}

const LABEL_PX = 11;
/** Each section's placed labels by its layout record (the composer's last placement), for hits on them (hit.ts). */
export const placedRingLabels = new WeakMap<object, RingLabels>();
// Enum-like members read once (a read per frame allocated under the test runner's module wrapper).
const INK = SkinToken.ink;
const QUIET = SkinToken.inkQuiet;
const STROKE = PathMode.stroke;
const MOVE = PathVerb.move;
const LINE = PathVerb.line;
const LEFT = TextAlign.left;
const RIGHT = TextAlign.right;
const MIDDLE = TextBaseline.middle;
/** Px between a leader's end and its text, and kept from the box's edge. */
const PAD = 3;
/** A label's numbers in `pts`: four leader points and the text's anchor. */
const STRIDE = 10;

/** The kept labels: per label its text as drawn, its side (1 right, -1 left), anchor and leader (CSS px). */
export class RingLabels {
  n = 0;
  text: string[] = [];
  side = new Int8Array(0);
  /**
   * Per label (`STRIDE` each): the leader's four points (on the item's ring, at the rings' edge on its angle, at the
   * labels' radius on the label's angle, at the text) then the text's anchor.
   */
  pts = new Float64Array(0);
  /** Labels left out (their side held no more lines, or no three graphemes fit): the composer reports them once. */
  hidden = 0;
  /** Per label its item and its track's point on the unrolled axis (`RingLabelInput`), and its text's width (CSS px). */
  item = new Int32Array(0);
  ux = new Float64Array(0);
  /** The layout version they were placed at (the composer's): a hit reads them only at that version. */
  at = -1;
  private readonly leaders = new PathBatch(32);
  private readonly right = new TextBatch(32);
  private readonly left = new TextBatch(32);

  /** Places `items` in `frame`, texts cut to the room their side leaves (edge-rate). */
  place(items: readonly RingLabelInput[], f: RingLabelFrame, text: TextManagerApi): this {
    const font = text.font('ui', LABEL_PX);
    const ellipsis = text.label('…', font).width;
    const kept: { i: RingLabelInput; side: number; y: number }[] = [];
    const cap = Math.max(0, Math.floor((f.bottom - f.top) / f.line));
    for (const side of [1, -1]) {
      const mine = items
        .filter((x) => (Math.cos(x.angle) >= 0 ? 1 : -1) === side)
        .sort((a, b) => a.priority - b.priority)
        .slice(0, cap)
        .map((i) => ({ i, side, y: f.cy + f.r * Math.sin(i.angle) }))
        .sort((a, b) => a.y - b.y);
      // A line apart: down from the top, then up from the bottom.
      for (let k = 0; k < mine.length; k++)
        mine[k].y = Math.max(mine[k].y, k > 0 ? mine[k - 1].y + f.line : f.top + f.line / 2);
      for (let k = mine.length - 1; k >= 0; k--)
        mine[k].y = Math.min(mine[k].y, k < mine.length - 1 ? mine[k + 1].y - f.line : f.bottom - f.line / 2);
      kept.push(...mine);
    }
    this.n = 0;
    this.hidden = items.length - kept.length;
    this.text = [];
    this.side = new Int8Array(kept.length);
    this.pts = new Float64Array(STRIDE * kept.length);
    this.item = new Int32Array(3 * kept.length);
    this.ux = new Float64Array(3 * kept.length);
    const fit = new Float64Array(2);
    fit[1] = ellipsis;
    for (const { i, side, y } of kept) {
      const dy = y - f.cy;
      // On the labels' circle at its y (beside the ring where the circle does not reach), within the box.
      const on = f.cx + side * Math.max(Math.sqrt(Math.max(0, f.r * f.r - dy * dy)), 12);
      const x = Math.min(f.right - 4 * PAD, Math.max(f.left + 4 * PAD, on));
      const p = text.label(i.text, font);
      fit[0] = side > 0 ? f.right - x - 2 * PAD : x - f.left - 2 * PAD;
      const m = fitIn(p, fit);
      if (m < Math.min(3, p.ends.length)) {
        this.hidden++;
        continue;
      }
      const k = this.n++;
      this.text.push(fitted(p, m));
      this.side[k] = side;
      this.item[3 * k] = i.layer ?? -1;
      this.item[3 * k + 1] = i.annotation ?? -1;
      this.item[3 * k + 2] = i.row ?? -1;
      this.ux[3 * k] = i.ux ?? 0;
      this.ux[3 * k + 1] = i.uy ?? 0;
      this.ux[3 * k + 2] = text.label(this.text[k], font).width;
      const c = Math.cos(i.angle);
      const s = Math.sin(i.angle);
      // The label's own angle, and a point on it short of the labels' circle, within the box.
      const la = Math.atan2(dy, x - f.cx);
      const lc = Math.cos(la);
      const ls = Math.sin(la);
      const room = ls < 0 ? (f.cy - f.top - PAD) / -ls : ls > 0 ? (f.bottom - PAD - f.cy) / ls : f.r;
      const elbow = Math.max(f.edge + PAD, Math.min(f.r - 2 * PAD, room, Math.hypot(x - f.cx, dy) - PAD));
      const o = STRIDE * k;
      this.pts.set(
        [
          f.cx + (i.r0 + 1) * c,
          f.cy + (i.r0 + 1) * s,
          f.cx + f.edge * c,
          f.cy + f.edge * s,
          f.cx + elbow * lc,
          f.cy + elbow * ls,
          x,
          y,
          x + side * PAD,
          y,
        ],
        o,
      );
    }
    return this;
  }

  /**
   * The label whose text holds surface point (x, y) (CSS px, a line high, `PAD` around it); -1: none. Read from `x`
   * and `y` of a record: a double passed to a call is boxed. Allocation-free.
   */
  labelAt(at: { readonly x: number; readonly y: number }): number {
    for (let k = 0; k < this.n; k++) {
      const o = STRIDE * k;
      const x0 = this.pts[o + 8];
      const w = this.ux[3 * k + 2] + PAD;
      const dx = this.side[k] > 0 ? at.x - x0 : x0 - at.x;
      if (dx >= -PAD && dx <= w && Math.abs(at.y - this.pts[o + 9]) <= LABEL_PX / 2 + 1) return k;
    }
    return -1;
  }

  /** The leaders (quiet ink, 1 px) and the texts (ink) into `be`, device px. Allocation-free. */
  draw(be: RenderBackend, text: TextManagerApi, dpr: number): void {
    if (this.n === 0) return;
    const leaders = this.leaders;
    leaders.reset();
    this.right.reset();
    this.left.reset();
    const font = text.font('ui', LABEL_PX);
    this.right.style(font, LEFT, MIDDLE);
    this.left.style(font, RIGHT, MIDDLE);
    this.right.begin(INK);
    this.left.begin(INK);
    // 1 CSS px at any DPR: the width is in half device px (VL-23).
    const w = Math.max(2, Math.round(2 * dpr)) | 0;
    for (let k = 0; k < this.n; k++) {
      const o = STRIDE * k;
      leaders.begin(QUIET, 255, STROKE, w);
      for (let j = 0; j < 4; j++) {
        const at = leaders.verb(j === 0 ? MOVE : LINE);
        leaders.coords[at] = this.pts[o + 2 * j] * dpr;
        leaders.coords[at + 1] = this.pts[o + 2 * j + 1] * dpr;
      }
      const t = this.side[k] > 0 ? this.right : this.left;
      const i = t.add(this.text[k]);
      t.x[i] = this.pts[o + 8] * dpr;
      t.y[i] = Math.round(this.pts[o + 9] * dpr);
    }
    be.paths(leaders);
    if (this.right.n > 0) be.text(this.right);
    if (this.left.n > 0) be.text(this.left);
  }
}
