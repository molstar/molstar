/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Anchor, NoColor, type PassWriter } from './draw-lists.js';
import { cellEdge } from './letters-draw.js';
import { drawColumnCells, drawTitles } from './rows-columns.js';
import { TogglePx } from './rows-header.js';
import type { Clades } from './rows-model.js';
import type { RowsLay } from './rows-lay.js';

// The rows layer's draw (_spec/README.md "Render", The rows layer): the cells, letters, pending rows and the pinned
// band's rule, the header (the tree's branches and wedges, labels, numbers, the chip, the row columns), all through the
// writer's `args` (nothing allocated once warm).

const AnchorStart: number = Anchor.start;
const AnchorEnd: number = Anchor.end;

// ---- drawing (nothing allocated once warm: numbers through the writer's `args`) ----

/** Clips to a box from its top `y` down to the band's foot (the scrolled rows below the pinned band). */
function clipBelow(lay: RowsLay, out: PassWriter, y: number) {
  const g = out.args;
  g[0] = lay.box.x;
  g[1] = y;
  g[2] = lay.box.width;
  g[3] = Math.max(0, lay.box.y + lay.box.height - y);
  out.clipArgs();
}

export function drawData(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  // the scrolled rows below the pinned band, clipped to it: the image tiles, the rects, their letters
  clipBelow(lay, out, lay.box.y + lay.pinnedHeight);
  if (lay.neutral) drawNeutral(lay, out);
  for (let k = 0; k < lay.images.length; k++) {
    g[0] = lay.ib[4 * k];
    g[1] = lay.ib[4 * k + 1];
    g[2] = lay.ib[4 * k + 2];
    g[3] = lay.ib[4 * k + 3];
    out.imageArgs(lay.images[k]);
  }
  drawRects(lay, out, 0, lay.pinnedFrom);
  if (lay.pending.length) drawPending(lay, out);
  if (lay.letters) drawLetters(lay, out, 0, lay.pinnedFrom);
  out.unclip();
  if (!lay.display.pinned) return;
  // the pinned rows on their own, over a band that does not scroll
  clipBelow(lay, out, lay.box.y);
  drawRects(lay, out, lay.pinnedFrom, lay.n);
  if (lay.letters) drawLetters(lay, out, lay.pinnedFrom, lay.n);
  drawPinnedRule(lay, out);
  out.unclip();
}

function drawRects(lay: RowsLay, out: PassWriter, from: number, to: number) {
  const g = out.args;
  // a letter cell's hairline (`seq-cell-edge`: dark only, none in light), as the letters layer's
  const edge = lay.letters ? lay.tokens.colors['seq-cell-edge'] : NoColor;
  const px = edge === NoColor ? 0 : 1 / lay.axis.dpr;
  for (let k = from; k < to; k++) {
    g[0] = lay.x0[k];
    g[1] = lay.y[k];
    g[2] = lay.x1[k] - lay.x0[k];
    g[3] = lay.h[k];
    g[4] = lay.fill[k];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
    if (edge !== NoColor) {
      g[5] = edge;
      g[6] = px;
      cellEdge(out);
    }
  }
}

function drawLetters(lay: RowsLay, out: PassWriter, from: number, to: number) {
  const g = out.args;
  let y = Number.NaN;
  let open = false;
  for (let k = from; k < to; k++) {
    const code = lay.code[k];
    if (!code) continue;
    const base = lay.y[k] + lay.h[k] / 2 + lay.fontSize * 0.35;
    if (base !== y) {
      if (open) out.endGlyphs();
      g[1] = base;
      out.beginGlyphsArgs(lay.font, NoColor);
      open = true;
      y = base;
    }
    g[0] = (lay.x0[k] + lay.x1[k]) / 2;
    g[4] = lay.ink[k];
    out.letterArgs(code);
  }
  if (open) out.endGlyphs();
}

function drawPending(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  for (let k = 0; k + 1 < lay.pending.length; k += 2) {
    g[0] = lay.box.x;
    g[1] = lay.pending[k];
    g[2] = lay.box.width;
    g[3] = lay.pending[k + 1];
    g[4] = lay.tokens.colors['seq-pending'];
    g[5] = NoColor;
    g[6] = 0;
    g[7] = 0;
    g[8] = 0;
    out.rectArgs();
  }
}

/** Lazy rows below a pixel: the band a neutral fill (never blank, never a pending stripe that never resolves). */
function drawNeutral(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  g[0] = lay.box.x;
  g[1] = lay.box.y + lay.pinnedHeight;
  g[2] = lay.box.width;
  g[3] = Math.max(0, lay.box.height - lay.pinnedHeight);
  g[4] = lay.tokens.colors['seq-missing'];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** A hairline under the pinned band. */
function drawPinnedRule(lay: RowsLay, out: PassWriter) {
  const g = out.args;
  g[0] = lay.box.x;
  g[1] = lay.box.y + lay.pinnedHeight - 0.5;
  g[2] = lay.box.width;
  g[3] = 1;
  g[4] = lay.tokens.colors['msui-line'];
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** The header: the tree's branches and wedges, the labels, the first numbers, the chip; the last numbers after the data. */
export function drawHeader(lay: RowsLay, out: PassWriter) {
  const { box, config, tokens } = lay;
  const g = out.args;
  g[0] = config.inset;
  g[1] = box.y;
  g[2] = Math.max(0, box.x - config.inset);
  g[3] = box.height;
  out.clipArgs();
  const ink = tokens.colors['msui-text-muted'];
  const t = lay.tree;
  branchWidth(lay, 1);
  for (let k = 0; k + 3 < t.length; k += 4) branchRect(lay, out, k, ink);
  if (lay.wedges.length) drawWedges(lay, out, ink, -1, false);
  if (lay.toggleX >= 0) drawToggle(lay, out, ink);
  const text = tokens.colors['msui-text'];
  // the names inline, not `row-names.ts`'s `drawNames`: shared with the matrix's, the call site saw two array kinds and
  // boxed a number a frame (16 B, the rows' allocation tier at GATE 1); the layout side stays shared (`addName`)
  for (let k = 0; k < lay.labelText.length; k++) {
    g[0] = lay.labelX;
    g[1] = lay.labelY[k];
    g[2] = lay.labelW[k];
    g[3] = lay.labelMax;
    g[4] = text;
    out.textArgs(lay.labelText[k], AnchorStart, lay.labelFont);
  }
  for (let k = 0; k < lay.numberText.length; k++) {
    g[0] = lay.numberX;
    g[1] = lay.numberY[k];
    g[2] = lay.numberW[k];
    g[3] = lay.numberW[k];
    g[4] = ink;
    out.textArgs(lay.numberText[k], AnchorEnd, config.itemFont);
  }
  if (lay.chip) {
    g[0] = lay.chipX;
    g[1] = lay.chipY;
    g[2] = lay.chipW;
    g[3] = lay.chipW;
    g[4] = tokens.colors['msui-accent'];
    out.textArgs(lay.chip, AnchorStart, config.itemFont);
  }
  if (lay.cols.n) drawColumnCells(lay, out, 0);
  out.unclip();
  if (lay.cols.n) {
    // the right columns, clipped to the band as the rows are
    g[0] = box.x + box.width;
    g[1] = box.y;
    g[2] = Math.max(0, lay.cols.rightEnd - g[0]);
    g[3] = box.height;
    out.clipArgs();
    drawColumnCells(lay, out, 1);
    out.unclip();
  }
  if (lay.cols.titleN || lay.cols.chip) drawTitles(lay, out);
  // the last numbers in the right gutter, clipped to the band as its rows are (a partly shown row ends at the foot)
  if (lay.lastText.length) {
    g[0] = box.x + box.width;
    g[1] = box.y;
    g[2] = lay.numberRightX - g[0] + lay.lastMaxW;
    g[3] = box.height;
    out.clipArgs();
  }
  for (let k = 0; k < lay.lastText.length; k++) {
    g[0] = lay.numberRightX;
    g[1] = lay.lastY[k];
    g[2] = lay.lastW[k];
    g[3] = lay.lastW[k];
    g[4] = ink;
    out.textArgs(lay.lastText[k], AnchorStart, config.itemFont);
  }
  if (lay.lastText.length) out.unclip();
}

/** Sets a branch's width `t` (px) in whole device pixels into `BranchW[0]`, one at least (_spec/README.md "Render", The rows layer, The tree's look). */
export function branchWidth(lay: RowsLay, t: number) {
  const dpr = lay.axis.dpr;
  BranchW[0] = Math.max(1, Math.round(t * dpr)) / dpr;
}

/**
 * The branches' width `branchWidth` set, read by `branchRect` and `drawWedges`; a wedge's corner radius (`softCorners`):
 * fractional numbers through an array (one passed to a call that is not inlined is boxed: the hover pass allocates
 * nothing).
 */
const BranchW = new Float64Array(2);

/**
 * Tree segment `k` (of `lay.tree`) as a rect of thickness `t` in `ink`, its edges on device pixels: across the segment
 * its centre's pixel run of the branch's width, along it from each end less half the width to the end plus half (so a
 * branch and its node's vertical meet square, and every branch draws one width).
 */
export function branchRect(lay: RowsLay, out: PassWriter, k: number, ink: number, tr: ArrayLike<number> = lay.tree) {
  const g = out.args;
  const dpr = lay.axis.dpr;
  const w = BranchW[0];
  const vertical = tr[k] === tr[k + 2];
  const across = vertical ? tr[k] : tr[k + 1];
  const a = vertical ? Math.min(tr[k + 1], tr[k + 3]) : Math.min(tr[k], tr[k + 2]);
  const b = vertical ? Math.max(tr[k + 1], tr[k + 3]) : Math.max(tr[k], tr[k + 2]);
  const c0 = Math.round((across - w / 2) * dpr) / dpr;
  const a0 = Math.round((a - w / 2) * dpr) / dpr;
  const b1 = Math.round((b - w / 2) * dpr) / dpr + w;
  g[0] = vertical ? c0 : a0;
  g[1] = vertical ? a0 : c0;
  g[2] = vertical ? w : b1 - a0;
  g[3] = vertical ? b1 - a0 : w;
  g[4] = ink;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The tree's toggle (a narrow frame's collapsed tree): a small tree in the branches' ink and width, a stem and three
 * branches, in the glyph's box; snapped to the device grid as the branches are.
 */
function drawToggle(lay: RowsLay, out: PassWriter, ink: number) {
  const x = lay.toggleX;
  const y = lay.toggleY;
  const s = TogglePx;
  const t = ToggleSegments;
  // the segments as a tree's (x0, y0, x1, y1), drawn by the branches' own rect
  t[0] = x + 2;
  t[1] = y + 2;
  t[2] = x + 2;
  t[3] = y + s - 2;
  t[4] = x + 2;
  t[5] = y + 2;
  t[6] = x + s - 2;
  t[7] = y + 2;
  t[8] = x + 2;
  t[9] = y + s / 2;
  t[10] = x + s - 4;
  t[11] = y + s / 2;
  t[12] = x + 2;
  t[13] = y + s - 2;
  t[14] = x + s - 2;
  t[15] = y + s - 2;
  // pressed (the tree is open): an underline, the toggle folding it back
  t[16] = x;
  t[17] = y + s + 1;
  t[18] = x + s;
  t[19] = y + s + 1;
  for (let k = 0; k < (lay.toggleOpen ? 20 : 16); k += 4) branchRect(lay, out, k, ink, t);
}

/** The toggle's segments, reused. */
const ToggleSegments = new Float64Array(20);

/** A wedge's corners and its outline with each corner softened (3 corners, 4 points each: two cuts, two between). */
const WedgeCorners = new Float64Array(6);
const Wedge = new Float64Array(24);

/**
 * The folded clades' wedges (_spec/README.md "Render", The rows layer, The tree's look): a triangle from the clade's
 * node to its deepest leaf's x, its far side half the row high, its corners softened by the shape rule, filled with a
 * light tint of `ink` (a step darker `hovered`) and outlined in `ink` at the branches' width `t`; all, or (`clade` ≥ 0)
 * those inside that clade.
 */
export function drawWedges(lay: RowsLay, out: PassWriter, ink: number, clade: number, hovered: boolean) {
  const w = lay.wedges;
  const g = out.args;
  const cl = lay.clades;
  const { config } = lay;
  const half = lay.rowScale / 4;
  const stroke = BranchW[0];
  // a tint of `ink` at a share of its own opacity
  const fill = ((ink & 0xffffff00) | Math.round((ink & 0xff) * (hovered ? 0.36 : 0.18))) >>> 0;
  for (let k = 0; k + 4 < w.length; k += 5) {
    if (clade >= 0 && cl && !within(cl, w[k + 4], clade)) continue;
    const x0 = w[k];
    const yc = (w[k + 1] + w[k + 3]) / 2;
    const x1 = w[k + 2];
    WedgeCorners[0] = x0;
    WedgeCorners[1] = yc;
    WedgeCorners[2] = x1;
    WedgeCorners[3] = yc - half;
    WedgeCorners[4] = x1;
    WedgeCorners[5] = yc + half;
    // the squircle rule's radius, inline (`squircleRadius` given fractional numbers boxes them)
    const len = x1 - x0;
    const r =
      len < config.squircleMin || 2 * half < config.squircleMin
        ? 0
        : Math.min(2 * half * config.squircleRatio, config.squircleMax, len / 2, half);
    BranchW[1] = r;
    const n = softCorners();
    g[4] = fill;
    g[5] = ink;
    g[6] = stroke;
    g[7] = 0;
    out.pathArgs(Wedge, n, true);
  }
}

/**
 * The triangle in `WedgeCorners` into `Wedge`, each corner cut `r` along its two sides (at most a third of the shorter)
 * and joined by a quadratic curve through the corner, sampled: its point count (3, the corners, when `r` is 0).
 */
function softCorners(): number {
  const r = BranchW[1];
  const p = WedgeCorners;
  if (!(r > 0)) {
    for (let i = 0; i < 6; i++) Wedge[i] = p[i];
    return 3;
  }
  let o = 0;
  for (let c = 0; c < 3; c++) {
    const px = p[2 * c];
    const py = p[2 * c + 1];
    const ax = p[2 * ((c + 2) % 3)];
    const ay = p[2 * ((c + 2) % 3) + 1];
    const bx = p[2 * ((c + 1) % 3)];
    const by = p[2 * ((c + 1) % 3) + 1];
    const la = Math.sqrt((ax - px) * (ax - px) + (ay - py) * (ay - py));
    const lb = Math.sqrt((bx - px) * (bx - px) + (by - py) * (by - py));
    const d = Math.min(r, la / 3, lb / 3);
    const sx = px + ((ax - px) * d) / la;
    const sy = py + ((ay - py) * d) / la;
    const ex = px + ((bx - px) * d) / lb;
    const ey = py + ((by - py) * d) / lb;
    for (let s = 0; s < 4; s++) {
      const u = s / 3;
      const v = 1 - u;
      Wedge[o++] = v * v * sx + 2 * u * v * px + u * u * ex;
      Wedge[o++] = v * v * sy + 2 * u * v * py + u * u * ey;
    }
  }
  return 12;
}

/** Node `k` is `c` or below it (its leaves are within `c`'s; a chain of one child shares them: by depth too). */
export function within(cl: Clades, k: number, c: number): boolean {
  return cl.lo[k] >= cl.lo[c] && cl.hi[k] <= cl.hi[c] && cl.depth[k] >= cl.depth[c];
}
