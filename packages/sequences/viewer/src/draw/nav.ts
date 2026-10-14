/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A navigator drawn from its parts (seqspec:viewer#extension-points), part of the renderer API: the whole axis on the
// nav pass (on a layout, size or skin change, never a pan: the map is the fit's), the window, the channels' marks and
// the grips on its overlay every frame. The axis spans exactly its blocks (`navX0s`, `navX1s`), its right edge the
// map's; marks per lane when stacked; on the ruler only `hairline` and `labels`. A miniature or overview is an image
// the composer keeps (`KeyedImage`). Device px as int32 throughout, batches owned and handed once.

import type { Axis, AxisView } from "../axis.ts";
import { MarkBit, type OverlayMarks } from "../marks.ts";
import { DefaultNavigatorParts, type NavigatorEntry, type NavigatorParts } from "../navigators.ts";
import { fitIn, fitted, type TextManagerApi } from "../text.ts";
import {
  type BackendImages,
  PathBatch,
  PathVerb,
  RectBatch,
  type RenderBackend,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from "./backend.ts";

/** What a navigator draws from: filled by the composer per pass (`from`, then the band and window), reused. */
export class NavView {
  backend: RenderBackend = null as unknown as RenderBackend;
  text: TextManagerApi = null as unknown as TextManagerApi;
  parts: NavigatorParts = DefaultNavigatorParts;
  dpr = 1;
  axis: Axis | null = null;
  /** The map's span from the data area's left (`navWidth`), and per block its extent (shared with the axis view). */
  left = 0;
  width = 0;
  x0s: Float64Array = new Float64Array(0);
  x1s: Float64Array = new Float64Array(0);
  /** The map's px per cell (`navCell`): zoom-independent. */
  cell = 1;
  /** The band: the navigator's row, or on the ruler the strip above its bottom border. CSS px. */
  top = 0;
  height = 0;
  /** The window as drawn and hit (windows.ts frame 8, 9); the pointer over the navigator (marks grow). */
  x0 = 0;
  x1 = 0;
  hot = false;
  /**
   * The window has nowhere to go (windows.ts `windowIdle`): drawn quiet, in its space, with no grips, ticks or hover.
   * A flag read per frame; the composer sets it from the viewport, so a wider view or more columns need no relayout.
   */
  idle = false;
  /** The grips lit (windows.ts `gripsLit`), and the band across they stay in. */
  lit = 0;
  gripTop = 0;
  gripBottom = 0;
  /** Per block its name (`labels`); what the channels light (marks and the hover's primary); null: none. */
  names: readonly string[] = [];
  marks: OverlayMarks | null = null;
  /** The axis in miniature or its overview, an image of the band's device size (`KeyedImage`); -1: none. */
  map = -1;
  /** The last axis drawn was the map: the window washes it outside and outlines itself in two tones. */
  overMap = false;
  /** The section's axis view under the current layout: the map, its blocks, the DPR. Allocation-free. */
  from(v: AxisView): this {
    this.axis = v.axis;
    this.dpr = v.dpr;
    this.left = v.left;
    this.width = v.navWidth;
    this.x0s = v.navX0s;
    this.x1s = v.navX1s;
    this.cell = v.navCell;
    return this;
  }
}

/** A navigator's draws: its own, else from its parts. */
export const navDraw = (e: NavigatorEntry, n: NavView) => (e.draw ? e.draw(n) : drawNavigator(n));
export const navDrawWindow = (e: NavigatorEntry, n: NavView) =>
  e.drawWindow ? e.drawWindow(n) : drawNavigatorWindow(n);

/** The window's wash over a map, as a byte, and its two-tone outline (light outside, dark inside). */
export const WINDOW_WASH = 158;
const WINDOW_LIGHT = 0xffffff;
/** A window's fill and edges (CG3: `--seq-window` and the grips' ink, never a selection's), read once. */
export const WINDOW_FILL = SkinToken.window;
export const WINDOW_EDGE_INK = SkinToken.grip;
const WINDOW_DARK = 0x111317;

const grounds = new PathBatch(16);
const under = new RectBatch(32);
const marks = new RectBatch(64);
const boxes = new RectBatch(16);
const labels = new TextBatch(16);
const outline = new RectBatch(8);
const grips = new PathBatch(8);
const img = new Float64Array(8);
/** The band in device px: y, height, lanes, a lane's height. */
const band = new Int32Array(4);
let shown = new Int32Array(32);
const fitArgs = new Float64Array(2);
/**
 * `drawGrips`' input: the window's ends along its axis (0, 1), the span the grips stay in (2, 3), their center across
 * (4) and the band across they stay in (5, 6), device px; the DPR (7). A double passed to a call is boxed: filled here.
 */
export const gripAt = new Float64Array(8);
/** `windowOverMap`'s input, device px: the window (x0, y0, x1, y1), then the map's box (x0, y0, x1, y1). */
export const windowBox = new Int32Array(8);
/** The grips' along and across lengths and least gap (CSS px). */
export const GRIP_ALONG = 3;
export const GRIP_ACROSS = 10;
const GRIP_GAP = 2;

/** The band (centered in the row; on the ruler at its bottom) and its lanes (a lane per block, stacked with `lanes`). */
function bandOf(n: NavView): void {
  const d = n.dpr;
  const th = Math.min(n.parts.thickness, n.height);
  const y = n.parts.onRuler ? n.top + n.height - th : n.top + (n.height - th) / 2;
  band[0] = Math.round(y * d) | 0;
  band[1] = Math.max(1, Math.round((y + th) * d) - band[0]) | 0;
  band[2] = n.parts.lanes && n.axis?.stacked ? n.axis.count : 1;
  band[3] = Math.max(1, (band[1] / band[2]) | 0);
}

/** The map image's size for this navigator, device px, into `out`: the data area's map across, the axis band down. */
export function mapSize(n: NavView, out: Int32Array): void {
  bandOf(n);
  out[0] = Math.max(1, Math.round((n.left + n.width) * n.dpr) - Math.round(n.left * n.dpr)) | 0;
  out[1] = band[1];
}

/** The whole axis: a ground per block (per lane), then the first of its `axis` parts that applies. */
export function drawNavigator(n: NavView): void {
  const p = n.parts;
  const ax = n.axis;
  n.overMap = false;
  if (!ax || n.backend === null) return;
  grounds.reset();
  under.reset();
  boxes.reset();
  labels.reset();
  const d = n.dpr;
  bandOf(n);
  const lanes = band[2];
  const lh = band[3];
  for (let b = 0; p.ground && band[1] > 2 && b < ax.count; b++) {
    const x = Math.round(n.x0s[b] * d) | 0;
    grounds.begin(SkinToken.surfaceAlt);
    const o = grounds.verb(PathVerb.roundRect);
    grounds.coords[o] = x;
    grounds.coords[o + 1] = band[0] + (lanes > 1 ? b * lh : 0);
    grounds.coords[o + 2] = Math.max(1, Math.round(n.x1s[b] * d) - x);
    grounds.coords[o + 3] = lanes > 1 ? lh : band[1];
    grounds.coords[o + 4] = Math.min(Math.round(p.radius * d), (lanes > 1 ? lh : band[1]) >> 1);
  }
  if (grounds.n > 0) n.backend.paths(grounds);
  for (let k = 0; k < p.axis.length; k++) {
    const a = p.axis[k];
    if (a === "overview" || a === "miniature") {
      if (p.onRuler || n.map < 0) continue;
      img[0] = img[1] = 0;
      img[2] = img[6] = Math.max(1, Math.round((n.left + n.width) * d) - Math.round(n.left * d));
      img[3] = img[7] = band[1];
      img[4] = Math.round(n.left * d);
      img[5] = band[0];
      n.backend.image(n.map, img, false);
      n.overMap = true;
      return;
    }
    if (a === "labels") {
      drawLabels(n);
      return;
    }
    if (a === "hairline") {
      const hair = Math.max(1, Math.round(d)) | 0;
      under.begin(SkinToken.inkFaint);
      for (let b = 0; b < ax.count; b++) {
        const x = Math.round(n.x0s[b] * d) | 0;
        const y = p.onRuler
          ? band[0] + band[1] - hair
          : lanes > 1
            ? band[0] + b * lh + ((lh - hair) >> 1)
            : band[0] + ((band[1] - hair) >> 1);
        under.push(x, y, Math.max(hair, Math.round(n.x1s[b] * d) - x), hair);
      }
      n.backend.rects(under);
      return;
    }
  }
}

/** Block names in outlined boxes (`band` with `axis: ["labels"]`). */
function drawLabels(n: NavView): void {
  const ax = n.axis as Axis;
  const tm = n.text;
  const d = n.dpr;
  const lanes = band[2];
  const lh = band[3];
  // 60% of the lane, 8 to 11 px: no Math builtin over a double (Maglev boxes its arguments: RN1).
  const px = ((lh * 3) / (5 * d)) | 0;
  const font = tm.font("ui", px < 8 ? 8 : px > 11 ? 11 : px);
  const lw = Math.max(1, Math.round(d)) | 0;
  const pad = Math.round(5 * d) | 0;
  const names = n.names;
  fitArgs[1] = tm.label("…", font).width;
  boxes.begin(SkinToken.inkFaint, 255, lw);
  labels.style(font, TextAlign.left, TextBaseline.middle);
  labels.begin(SkinToken.inkQuiet);
  for (let b = 0; b < ax.count; b++) {
    const x = Math.round(n.x0s[b] * d) | 0;
    const w = Math.max(2 * lw, Math.round(n.x1s[b] * d) - x) | 0;
    const y = band[0] + (lanes > 1 ? b * lh : 0);
    boxes.push(x, y, w, lanes > 1 ? lh : band[1]);
    const label = tm.label(b < names.length ? names[b] : "", font);
    fitArgs[0] = w / d - 10;
    const count = fitIn(label, fitArgs);
    if (count <= 0) continue;
    const i = labels.add(count >= label.ends.length ? label.text : fitted(label, count));
    labels.x[i] = x + pad;
    labels.y[i] = y + ((lanes > 1 ? lh : band[1]) >> 1);
  }
  n.backend.rects(boxes);
  if (labels.n > 0) n.backend.text(labels);
}

/**
 * The window: the outside faded (over a map: washed, the window outlined in two tones), an outlined box, a segment on
 * its light track (the rail's two tones), end ticks; the channels' marks per lane, the hover's primary; the end grips.
 */
export function drawNavigatorWindow(n: NavView): void {
  const p = n.parts;
  const ax = n.axis;
  if (!ax || n.backend === null) return;
  under.reset();
  marks.reset();
  const d = n.dpr;
  bandOf(n);
  const by = band[0];
  const bh = band[1];
  const x0 = Math.round(n.x0 * d) | 0;
  const x1 = Math.max(x0 + 1, Math.round(n.x1 * d)) | 0;
  const lo = Math.round(n.left * d) | 0;
  // The map's right edge, not the last block's (stacked, the widest lane's).
  const hi = Math.round((n.left + n.width) * d) | 0;
  const rowTop = Math.round(n.top * d) | 0;
  const rowBottom = Math.round((n.top + n.height) * d) | 0;
  const lw = Math.max(1, Math.round(d)) | 0;
  // Idle (nowhere to go): the window is the whole axis, so no wash, box or outline.
  if (!n.idle && n.overMap && (p.box || p.fade > 0)) {
    windowBox[0] = x0;
    windowBox[1] = rowTop;
    windowBox[2] = x1;
    windowBox[3] = rowBottom;
    windowBox[4] = lo;
    windowBox[5] = rowTop;
    windowBox[6] = hi;
    windowBox[7] = rowBottom;
    windowOverMap(n.backend, lw, SkinToken.surface, p.fade > 0 ? Math.round(p.fade * 255) : WINDOW_WASH);
  } else if (!n.idle && p.fade > 0) {
    under.begin(SkinToken.surface, Math.round(p.fade * 255));
    if (x0 > lo) under.push(lo, by, x0 - lo, bh);
    if (hi > x1) under.push(x1, by, hi - x1, bh);
  }
  if (p.box && !n.overMap && !n.idle && Math.min(x1, hi) > Math.max(x0, lo)) {
    const bx = Math.max(x0, lo);
    under.begin(WINDOW_FILL);
    under.push(bx, rowTop, Math.min(x1, hi) - bx, rowBottom - rowTop);
    under.begin(WINDOW_EDGE_INK, 255, lw);
    under.push(bx, rowTop, Math.min(x1, hi) - bx, rowBottom - rowTop);
  }
  // Across the band: the segment's center, else the band's.
  let mid = by + (bh >> 1);
  if (p.segment > 0) {
    const h = Math.max(1, Math.round((n.hot && !n.idle && p.hoverGrow ? p.segment + 2 : p.segment) * d)) | 0;
    const y = p.onRuler ? by + bh - h : by + ((bh - h) >> 1);
    mid = y + (h >> 1);
    under.begin(SkinToken.inkFaint);
    for (let b = 0; b < ax.count; b++) {
      const x = Math.round(n.x0s[b] * d) | 0;
      under.push(x, y, Math.max(1, Math.round(n.x1s[b] * d) - x), h);
    }
    // Idle: the window is the quiet track itself, in the disabled tone (`inkFaint`), so there is no second line.
    if (!n.idle) {
      under.begin(SkinToken.inkQuiet);
      under.push(x0, y, x1 - x0, h);
    }
  }
  if (p.ticks && !p.box && !n.idle) {
    under.begin(WINDOW_EDGE_INK);
    under.push(x0, by, lw, bh);
    under.push(x1 - lw, by, lw, bh);
  }
  if (under.n > 0) n.backend.rects(under);
  if (n.marks) {
    if (p.highlightMark) {
      marksOf(n, MarkBit.highlight, SkinToken.highlight);
      marksOf(n, MarkBit.accent, SkinToken.accent);
      marksOf(n, MarkBit.selection, SkinToken.selectionBorder);
    }
    if (p.hoverMark) primaryOf(n);
    if (marks.n > 0) n.backend.rects(marks);
  }
  if (n.idle) return;
  gripAt[0] = x0;
  gripAt[1] = x1;
  gripAt[2] = lo;
  gripAt[3] = hi;
  gripAt[4] = mid;
  gripAt[5] = Math.round(n.gripTop * d) | 0;
  gripAt[6] = Math.round(n.gripBottom * d) | 0;
  gripAt[7] = d;
  drawGrips(n.backend, false, n.lit);
}

/** One channel's runs as marks along the bottom of each block's lane (2 px), through the axis's run converter. */
function marksOf(n: NavView, bit: number, color: number): void {
  const ax = n.axis as Axis;
  const m = (n.marks as OverlayMarks).runs(bit);
  const d = n.dpr;
  const lanes = band[2];
  const mh = Math.max(1, Math.round(2 * d)) | 0;
  marks.begin(color);
  for (let b = 0; b < ax.count; b++) {
    const to = m.to(b);
    if (to <= m.from(b) || ax.shown[b] === 0) continue;
    const per = (n.x1s[b] - n.x0s[b]) / ax.shown[b];
    const y = (lanes > 1 ? band[0] + (b + 1) * band[3] : band[0] + band[1]) - mh;
    for (let r = m.from(b); r < to; r++) {
      let c = ax.runToShown(b, m.start[r], m.stop[r], shown);
      if (2 * c > shown.length) {
        shown = new Int32Array(4 * c);
        c = ax.runToShown(b, m.start[r], m.stop[r], shown);
      }
      for (let q = 0; q < c; q++) {
        const xa = Math.round((n.x0s[b] + shown[2 * q] * per) * d) | 0;
        const xb = Math.round((n.x0s[b] + shown[2 * q + 1] * per) * d) | 0;
        marks.push(xa, y, Math.max(1, xb - xa), mh);
      }
    }
  }
}
/** The hover's primary as a line across its block's lane: at the center of its first lit position the layout shows. */
function primaryOf(n: NavView): void {
  const ax = n.axis as Axis;
  const mk = n.marks as OverlayMarks;
  const d = n.dpr;
  const lanes = band[2];
  const w = Math.max(1, Math.round(d)) | 0;
  marks.begin(SkinToken.accent);
  for (let b = 0; b < ax.count; b++) {
    const i = mk.primary(b);
    const k = i < 0 || ax.shown[b] === 0 ? -1 : ax.shownOf(b, i);
    if (k < 0) continue;
    const per = (n.x1s[b] - n.x0s[b]) / ax.shown[b];
    const x = Math.round((n.x0s[b] + (k + 0.5) * per) * d) | 0;
    marks.push(x, lanes > 1 ? band[0] + b * band[3] : band[0], w, lanes > 1 ? band[3] : band[1]);
  }
}

/**
 * Over a map (a navigator's miniature, the vertical overview): the map outside the window washed with `ground` at
 * `alpha` (a byte), the window at full color outlined in two tones, a light line `lw` device px wide outside a dark
 * one, inside the map's box. A window thinner than both lines on either side is outlined at their width, centered on
 * it. From `windowBox` (device px).
 */
export function windowOverMap(backend: RenderBackend, lw: number, ground: number, alpha: number): void {
  const w = windowBox;
  const bx0 = w[4];
  const by0 = w[5];
  const bx1 = w[6];
  const by1 = w[7];
  const least = 4 * lw;
  if (bx1 - bx0 < least || by1 - by0 < least) return;
  let x0 = Math.max(w[0], bx0);
  let x1 = Math.min(w[2], bx1);
  let y0 = Math.max(w[1], by0);
  let y1 = Math.min(w[3], by1);
  if (x1 - x0 < least) {
    x0 = Math.max(bx0, Math.min(bx1 - least, (x0 + x1 - least) >> 1));
    x1 = x0 + least;
  }
  if (y1 - y0 < least) {
    y0 = Math.max(by0, Math.min(by1 - least, (y0 + y1 - least) >> 1));
    y1 = y0 + least;
  }
  outline.reset();
  if (alpha > 0) {
    outline.begin(ground, alpha);
    if (y0 > by0) outline.push(bx0, by0, bx1 - bx0, y0 - by0);
    if (by1 > y1) outline.push(bx0, y1, bx1 - bx0, by1 - y1);
    if (x0 > bx0) outline.push(bx0, y0, x0 - bx0, y1 - y0);
    if (bx1 > x1) outline.push(x1, y0, bx1 - x1, y1 - y0);
  }
  outline.begin(WINDOW_LIGHT, 242, lw);
  outline.push(x0, y0, x1 - x0, y1 - y0);
  outline.begin(WINDOW_DARK, 217, lw);
  outline.push(x0 + lw, y0 + lw, x1 - x0 - 2 * lw, y1 - y0 - 2 * lw);
  backend.rects(outline);
}

/** The pills in device px: start and end along, across, along and across lengths, radius, halo, band (7, 8), span (9, 10). */
const gp = new Int32Array(11);
/** Grip `which` (0 the start's, 1 the end's) as a rounded rect, grown by `grow` and kept inside the span and band. */
function pill(vertical: boolean, which: number, grow: number): void {
  const a = Math.max(gp[9], gp[which] - grow);
  const al = Math.min(gp[10], gp[which] + gp[3] + grow) - a;
  const c = Math.max(gp[7], gp[2] - grow);
  const cl = Math.min(gp[8], gp[2] + gp[4] + grow) - c;
  const o = grips.verb(PathVerb.roundRect);
  grips.coords[o] = vertical ? c : a;
  grips.coords[o + 1] = vertical ? a : c;
  grips.coords[o + 2] = vertical ? cl : al;
  grips.coords[o + 3] = vertical ? al : cl;
  grips.coords[o + 4] = Math.min(gp[5] + grow, al >> 1, cl >> 1);
}
/**
 * A window's two end grips, from `gripAt`: pills `GRIP_ALONG` along the window's axis (y with `vertical`) and
 * `GRIP_ACROSS` across, centered on its ends and kept inside the span and band, so at the full view they sit at the
 * span's far ends; a window shorter than both and their gap draws them apart, centered on it. Each on a
 * halo of the ground; `lit` bits (windows.ts `gripsLit`: 1 the start's, 2 the end's) in `gripHover`, the rest `grip`.
 */
export function drawGrips(backend: RenderBackend, vertical: boolean, lit: number): void {
  const a = gripAt;
  const d = a[7];
  const t = Math.max(2, Math.round(GRIP_ALONG * d)) | 0;
  const gap = Math.round(GRIP_GAP * d) | 0;
  const lo = a[2] | 0;
  const hi = a[3] | 0;
  const c0 = a[5] | 0;
  const c1 = a[6] | 0;
  const len = Math.min(Math.round(GRIP_ACROSS * d), c1 - c0) | 0;
  if (len < 2 || hi - lo < 2 * t + gap) return;
  let s = (a[0] | 0) - (t >> 1);
  let e = (a[1] | 0) - (t >> 1);
  if (e - s < t + gap) {
    s = ((s + e) >> 1) - ((t + gap) >> 1);
    e = s + t + gap;
  }
  s = Math.max(lo, Math.min(hi - 2 * t - gap, s));
  e = Math.max(s + t + gap, Math.min(hi - t, e));
  gp[0] = s;
  gp[1] = e;
  gp[2] = Math.max(c0, Math.min(c1 - len, (a[4] | 0) - (len >> 1)));
  gp[3] = t;
  gp[4] = len;
  gp[5] = t >> 1;
  gp[6] = Math.max(1, Math.round(d));
  gp[7] = c0;
  gp[8] = c1;
  gp[9] = lo;
  gp[10] = hi;
  grips.reset();
  for (let k = 0; k < 2; k++) {
    grips.begin(SkinToken.surface);
    pill(vertical, k, gp[6]);
  }
  for (let k = 0; k < 2; k++) {
    grips.begin((lit >> k) & 1 ? SkinToken.gripHover : SkinToken.grip);
    pill(vertical, k, 0);
  }
  backend.paths(grips);
}

/**
 * An image built at edge rate and kept while its key holds: the `seqv:map` miniature, a navigator's map, the vertical
 * overview's. The caller writes the key's numbers into `next` (a layout version, a size, the DPR, a skin stamp, data
 * versions), asks `stale`, and on true builds the image and `keep`s it (-1: nothing to draw, kept as such, so a frame
 * with nothing to draw builds nothing either); `builds` counts the images made.
 */
export class KeyedImage {
  id = -1;
  builds = 0;
  readonly next: Float64Array;
  private readonly key: Float64Array;
  private images: BackendImages | null = null;
  constructor(size = 8) {
    this.next = new Float64Array(size);
    this.key = new Float64Array(size).fill(Number.NaN);
  }
  /** True when nothing is kept for these images, or `next` differs from the key it was kept under. */
  stale(images: BackendImages): boolean {
    if (images !== this.images) return true;
    for (let k = 0; k < this.key.length; k++) if (this.key[k] !== this.next[k]) return true;
    return false;
  }
  /** A new image (or -1: none) under `next`; the previous one disposed. */
  keep(images: BackendImages, id: number): void {
    this.dispose();
    this.images = images;
    this.id = id;
    this.key.set(this.next);
    if (id >= 0) this.builds++;
  }
  dispose(): void {
    if (this.id >= 0) this.images?.dispose(this.id);
    this.id = -1;
    this.images = null;
    this.key.fill(Number.NaN);
  }
}
