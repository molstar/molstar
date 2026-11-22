/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Location, LocationPart } from '@molstar/seqx-core/query';
import type { AxisTransform } from './axis.js';
import { NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { walkStripes } from './highlight-stripes.js';
import { type Emphasis, type EmphasisHighlight } from './layer-type.js';
import { partOn, rectOver } from './marks.js';
import { CrossEdge, CrossTint, type GridLay } from './grid.js';

// The grid layer's emphasis (_spec/README.md "Render", The grid layer): the y path's rows at each lit y, the
// crosshair, a selected rect, the cursor's cell; the strip paints x.

/**
 * What the emphasis paints with: the fill, the stroke and its width (`Ink[0..2]`), set before each paint. The colors
 * (above 2^30, no small integers) and the widths go through a typed array: a number passed to a call that is not
 * inlined is boxed, and an emphasis redraw allocates nothing once warm.
 */
export const Ink = new Float64Array(3);

/**
 * Sets `Ink` from `Src` (a fill color and its opacity, a stroke color and its opacity, a width): a color at an opacity
 * (0 to 1), none at 0. The colors are read from an array, never passed (a color above 2^30 passed to a call is boxed).
 */
const Src = new Float64Array(5);

function ink() {
  const f = Src[0];
  const s = Src[2];
  Ink[0] = Src[1] > 0 ? ((f & 0xffffff00) | Math.round(Math.min(1, Src[1]) * (f & 0xff))) >>> 0 : NoColor;
  Ink[1] = Src[3] > 0 ? ((s & 0xffffff00) | Math.round(Math.min(1, Src[3]) * (s & 0xff))) >>> 0 : NoColor;
  Ink[2] = Src[4];
}

/** A range's extent on an axis, its first cell to its last, into `Box[o]` (low) and `Box[o + 1]` (high); false: hidden. */
const Box4 = new Float64Array(4);

function span(a: AxisTransform, s: number, e: number, o: number): boolean {
  if (!a.cell(s)) return false;
  let lo = a.sx0 < a.sx1 ? a.sx0 : a.sx1;
  let hi = a.sx0 < a.sx1 ? a.sx1 : a.sx0;
  if (!a.cell(e - 1)) return false;
  if (a.sx0 < lo) lo = a.sx0;
  if (a.sx1 < lo) lo = a.sx1;
  if (a.sx0 > hi) hi = a.sx0;
  if (a.sx1 > hi) hi = a.sx1;
  Box4[o] = lo;
  Box4[o + 1] = hi;
  return true;
}

/** A rect's x, y, width and height, set before `rect`. */
export const R = new Float64Array(4);

/** The rect's squircle radius and round ends (`rounded` sets them; 0: square, the crosshair's rows and columns). */
const RR = new Float64Array(2);

/** The rect `R` holds in `Ink`, its corners `RR`'s. */
export function rect(out: PassWriter) {
  const g = out.args;
  g[0] = R[0];
  g[1] = R[1];
  g[2] = R[2];
  g[3] = R[3];
  g[4] = Ink[0];
  g[5] = Ink[1];
  g[6] = Ink[2];
  g[7] = RR[0];
  g[8] = RR[1];
  out.rectArgs();
}

/**
 * `R` as a squircle into `RR` (every emphasis a squircle: a cell's outline, a selected rect): the squircle rule over
 * its size, square below the minimum, a side past the plot (x's data box, the band's top or bottom) cut square.
 */
function rounded(lay: GridLay) {
  const c = lay.config!;
  const a = lay.axis!;
  const b = lay.box;
  const w = R[2];
  const h = R[3];
  RR[0] = 0;
  RR[1] = 0;
  if (w < c.squircleMin || h < c.squircleMin) return;
  let round = (R[0] >= a.x0 - 0.01 ? 1 : 0) | (R[0] + w <= a.x0 + a.width + 0.01 ? 2 : 0);
  if (R[1] < b.y - 0.01) round |= 4;
  if (R[1] + h > b.y + b.height + 0.01) round |= 8;
  if ((round & 3) === 0) return;
  RR[0] = Math.min(h * c.squircleRatio, c.squircleMax, w / 2, h / 2);
  RR[1] = round;
}

/**
 * Rows across the plot at a part's ranges on y (each range's span down the plot), in `Ink`: the y half of the
 * crosshair, a selection's, a highlight's (_spec/README.md "Render", The grid layer).
 */
function paintRows(lay: GridLay, out: PassWriter, ranges: Int32Array) {
  const ay = lay.axisY!;
  const a = lay.axis!;
  const b = lay.box;
  // a stripe across the plot (a highlight's, where channels overlap): its edges on device pixels, shared with its
  // neighbors; a band is a bar along its leading edge
  stripeEdges(a);
  const w = RowX[1] - RowX[0];
  if (w <= 0) return;
  for (let r = 0; r + 1 < ranges.length; r += 2) {
    if (!span(ay, ranges[r], ranges[r + 1], 2)) continue;
    const top = Math.max(b.y, Box4[2]);
    const bottom = Math.min(b.y + b.height, Box4[3]);
    if (bottom <= top) continue;
    // across the plot: x's extent (the plot may be narrower than its band)
    R[0] = RowX[0];
    R[1] = top;
    R[2] = RowX[1] - RowX[0];
    R[3] = bottom - top;
    rect(out);
  }
}

/**
 * The x edges `[RowX[0], RowX[1]]` of stripe `slot` of `count` across the plot, on device pixels (shared with its
 * neighbors); a bar (`RowLook[2]` px) is along the leading edge. Through an array: a double returned or kept is boxed.
 */
const RowX = new Float64Array(2);

function stripeEdges(a: AxisTransform) {
  const slot = RowLook[0];
  const count = RowLook[1];
  const dpr = a.dpr;
  if (count > 1) {
    RowX[0] = Math.round((a.x0 + (slot * a.width) / count) * dpr) / dpr;
    RowX[1] = Math.round((a.x0 + ((slot + 1) * a.width) / count) * dpr) / dpr;
  } else {
    RowX[0] = a.x0;
    RowX[1] = a.x0 + a.width;
  }
  if (RowLook[2] > 0)
    RowX[1] = RowX[0] + Math.min(Math.max(Math.round(RowLook[2] * dpr) / dpr, 1 / dpr), RowX[1] - RowX[0]);
}

/** A stripe's look for `paintRows`: its slot, the stripe count, and a bar's width (0: the whole stripe). */
const RowLook = new Float64Array([0, 1, 0]);

/** The rows' highlight channels (`walkStripes`' context: set before it, since it is given no closure). */
const GParts: (LocationPart | undefined)[] = [];
const GIndex: number[] = [];
const GOutlines: number[] = [];
const GRun = new Int32Array(2);
let GLay!: GridLay;
let GOut!: PassWriter;
let GList!: readonly EmphasisHighlight[];

/**
 * The highlight channels' rows (the plot's y part; the x part is the strip's, drawn with the tracks'): where channels
 * overlap the row is split into equal stripes across the plot, in channel order (as in tracks, a band's stripes
 * perpendicular to its axis); `band` is a bar along a stripe's left edge; `outline` takes no stripe, one outline over
 * the whole row.
 */
function drawRowHighlights(lay: GridLay, out: PassWriter, hs: readonly EmphasisHighlight[]) {
  const config = lay.config!;
  let n = 0;
  let o = 0;
  for (let i = 0; i < hs.length; i++) {
    const part = partOn(hs[i].location, lay.ySystem);
    if (!part) continue;
    if (hs[i].style === 'outline') GOutlines[o++] = i;
    else {
      GParts[n] = part;
      GIndex[n] = i;
      n++;
    }
  }
  GLay = lay;
  GOut = out;
  GList = hs;
  if (n === 1) {
    rowsOf(lay, out, hs[GIndex[0]], GParts[0]!.ranges, 0, 1);
  } else if (n > 1) walkStripes(GParts, n, Math.max(1, config.emphasis.highlightStripes | 0), emitRows);
  for (let k = 0; k < o; k++) {
    const h = hs[GOutlines[k]];
    rowsOf(lay, out, h, partOn(h.location, lay.ySystem)!.ranges, 0, 1);
  }
  for (let c = 0; c < n; c++) GParts[c] = undefined;
}

function emitRows(channel: number, start: number, end: number, slot: number, count: number, _first: number) {
  GRun[0] = start;
  GRun[1] = end;
  rowsOf(GLay, GOut, GList[GIndex[channel]], GRun, slot, count);
}

/** A channel's ranges as rows in its style, in stripe `slot` of `count`. */
function rowsOf(lay: GridLay, out: PassWriter, h: EmphasisHighlight, ranges: Int32Array, slot: number, count: number) {
  const em = lay.config!.emphasis;
  const style = h.style ?? 'tint';
  Src[0] = h.color;
  Src[1] = style === 'tint' ? em.highlightTint : style === 'band' ? 1 : 0;
  Src[2] = h.color;
  Src[3] = style === 'outline' ? 1 : 0;
  Src[4] = style === 'outline' ? lay.config!.outline.selection : 0;
  ink();
  RowLook[0] = slot;
  RowLook[1] = count;
  RowLook[2] = style === 'band' ? em.highlightBand : 0;
  paintRows(lay, out, ranges);
  RowLook[0] = 0;
  RowLook[1] = 1;
  RowLook[2] = 0;
}

/**
 * A cell's box (x position `Pos[0]`, y position `Pos[1]`) in `Ink`; nothing when it is not shown. The positions through
 * an array: a hover's may be held as doubles, and a double passed to a call is boxed.
 */
const Pos = new Float64Array(2);

function paintCell(lay: GridLay, out: PassWriter) {
  const a = lay.axis!;
  const ay = lay.axisY!;
  const i = Pos[0];
  const j = Pos[1];
  if (!a.cell(i)) return;
  R[0] = Math.min(a.sx0, a.sx1);
  R[2] = Math.abs(a.sx1 - a.sx0);
  if (!ay.cell(j)) return;
  R[1] = Math.min(ay.sx0, ay.sx1);
  R[3] = Math.abs(ay.sx1 - ay.sx0);
  rounded(lay);
  rect(out);
  RR[0] = 0;
  RR[1] = 0;
}

/** A column down the plot at x position `Pos[0]` in `Ink` (the crosshair's x half, on the hover's own section). */
function paintColumn(lay: GridLay, out: PassWriter) {
  const a = lay.axis!;
  const i = Pos[0];
  if (!a.cell(i)) return;
  const b = lay.box;
  R[0] = Math.min(a.sx0, a.sx1);
  R[1] = b.y;
  R[2] = Math.abs(a.sx1 - a.sx0);
  R[3] = b.height;
  rect(out);
}

/** Each rect over the plot's two systems in `Ink` (a selected rectangle). */
function paintRects(lay: GridLay, out: PassWriter, loc: Location) {
  const a = lay.axis!;
  const ay = lay.axisY!;
  for (let k = 0; k < loc.rects.length; k++) {
    const r = loc.rects[k];
    if (r.a.system !== lay.system || r.b.system !== lay.ySystem) continue;
    for (let p = 0; p + 1 < r.a.ranges.length; p += 2)
      for (let q = 0; q + 1 < r.b.ranges.length; q += 2) {
        if (!span(a, r.a.ranges[p], r.a.ranges[p + 1], 0) || !span(ay, r.b.ranges[q], r.b.ranges[q + 1], 2)) continue;
        R[0] = Box4[0];
        R[1] = Box4[2];
        R[2] = Box4[1] - Box4[0];
        R[3] = Box4[3] - Box4[2];
        rounded(lay);
        rect(out);
        RR[0] = 0;
        RR[1] = 0;
      }
  }
}

/**
 * The emphasis on y (the strip paints x): highlights' rows; the selection's rows and rects (its own, else projected,
 * lighter); the hover's crosshair on this section (a column and a row and the cell's outline) or a projected hover's
 * rows; the cursor's cell when it is on a cell of the plot.
 */
export function drawEmphasis(lay: GridLay, pass: PassName, out: PassWriter, e: Emphasis) {
  const tokens = lay.tokens!;
  const config = lay.config!;
  const em = config.emphasis;
  switch (pass) {
    case 'highlight': {
      if (e.highlights && em.selectionPlace === 'over') drawRowHighlights(lay, out, e.highlights);
      return;
    }
    case 'selection': {
      const accent = tokens.colors['msui-accent'];
      const own = partOn(e.selection, lay.ySystem);
      const part =
        own ??
        (rectOver(e.selection, lay.system, lay.ySystem) ? undefined : partOn(e.projected?.selection, lay.ySystem));
      const k = own ? 1 : em.projected;
      if (part) {
        Src[0] = accent;
        Src[1] = em.selectionPlace === 'over' ? em.selectionTint * k : 0;
        Src[2] = accent;
        Src[3] = em.selectionOutline * k;
        Src[4] = config.outline.selection;
        ink();
        paintRows(lay, out, part.ranges);
      }
      if (e.selection?.rects.length) {
        Src[0] = accent;
        Src[1] = em.selectionTint;
        Src[2] = accent;
        Src[3] = em.selectionOutline;
        Src[4] = config.outline.selection;
        ink();
        paintRects(lay, out, e.selection);
      }
      return;
    }
    case 'hover': {
      const hover = tokens.colors['seq-hover'];
      const h = e.hover;
      if (
        h &&
        h.section === lay.section &&
        h.yPosition !== undefined &&
        h.position !== undefined &&
        h.element === undefined &&
        (h.ySystem === lay.ySystem || h.system === lay.system)
      ) {
        // on a grid over blocks the crosshair crosses every strip of the hovered x block (its column) and of its y
        // block (its row), the cell outlined on the one strip of both; the column and row tinted and edged in the
        // hover's ink (its contrast token: dark in light, light in dark), so they read on either end of a palette
        Src[0] = hover;
        Src[1] = CrossTint;
        Src[2] = hover;
        Src[3] = CrossEdge;
        Src[4] = config.outline.selection;
        ink();
        Pos[0] = h.position;
        Pos[1] = h.yPosition;
        if (h.system === lay.system) paintColumn(lay, out);
        Cell[0] = h.yPosition;
        Cell[1] = h.yPosition + 1;
        if (h.ySystem === lay.ySystem) paintRows(lay, out, Cell);
        if (h.ySystem !== lay.ySystem || h.system !== lay.system) return;
        Src[0] = 0;
        Src[1] = 0;
        Src[2] = hover;
        Src[3] = 1;
        Src[4] = config.outline.hover;
        ink();
        paintCell(lay, out);
        return;
      }
      const p = e.projected;
      if (p?.hover && p.hoverSection !== lay.sectionAddress) {
        const part = partOn(p.hover, lay.ySystem);
        if (!part) return;
        Src[0] = 0;
        Src[1] = 0;
        Src[2] = hover;
        Src[3] = em.projected;
        Src[4] = config.outline.hover;
        ink();
        paintRows(lay, out, part.ranges);
      }
      return;
    }
    case 'cursor': {
      const c = e.cursor;
      if (!c || c.ySystem !== lay.ySystem || c.yPosition === undefined || c.system !== lay.system) return;
      // on its own section (another of these systems marks its column)
      if (c.section !== undefined && c.section !== lay.sectionAddress) return;
      Src[0] = 0;
      Src[1] = 0;
      Src[2] = tokens.colors[em.cursorInk];
      Src[3] = 1;
      Src[4] = config.outline.cursor;
      ink();
      Pos[0] = c.position;
      Pos[1] = c.yPosition;
      paintCell(lay, out);
      return;
    }
  }
}

/** One cell's y range, reused (a hover's row). */
const Cell = new Int32Array(2);
