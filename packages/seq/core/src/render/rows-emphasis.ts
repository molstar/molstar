/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { NoColor, type PassWriter, Round } from './draw-lists.js';
import type { Emphasis } from './layer-type.js';
import { hoverColumnCells } from './rows-columns.js';
import { squircleRadius } from './squircle.js';
import type { RowsLay } from './rows-lay.js';
import { branchRect, branchWidth, drawWedges, within } from './rows-draw.js';

// The rows layer's emphasis along y (_spec/README.md "Render", The rows layer, Emphasis along y; Hover): the hovered
// row, cell and clade, the selection's rows and cells, the keyboard cursor, a row column's cell read whole (`tip`).

/**
 * A row column's cell cut short, read whole in the `tip` pass (_spec/README.md "Render", The rows layer, Emphasis along
 * y): the one cell under the pointer (a header hit naming it), never from a hover on the row's data or label; the
 * keyboard cursor on a row's header reads all of its row's. A function of its own (measured: no boxing).
 */
export function hoverColumns(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  if (!lay.cols.n) return;
  const h = emphasis?.hover;
  if (h && h.layer === lay.key && h.target === 'row' && h.state !== undefined && h.position !== undefined) {
    // the copy the pointer is on (a pinned row kept in place is shown twice): the last hit's, when it shows that row
    const r = lay.hitRow;
    const i =
      r >= 0 && r < lay.display.count && lay.display.rowAt(r) === h.position ? r : lay.display.displayRow(h.position);
    if (i >= 0) hoverColumnCells(lay, out, i, h.state);
  }
  const c = emphasis?.cursor;
  if (!c || c.part !== 'header' || c.ySystem !== lay.rowsSystem || c.yPosition === undefined) return;
  // the layer's address is its section's and more: `startsWith` and a code read allocate nothing
  const s = c.section;
  if (s !== undefined && !(lay.key.startsWith(s) && lay.key.charCodeAt(s.length) === 47)) return;
  // on every copy the cursor is drawn on (as `drawCursor` outlines both of a pinned row's)
  const i = lay.cursorRow(c.yPosition);
  if (i >= 0) hoverColumnCells(lay, out, i, undefined);
  const own = lay.display.placeRow(c.yPosition);
  if (own >= 0 && own !== i) hoverColumnCells(lay, out, own, undefined);
}

/**
 * A light band over `rows` display rows from `i` (one by default) as far as it is in the band: the scrolled rows
 * (`i` past the pinned band) never cover the pinned band.
 */
function hoverBand(lay: RowsLay, out: PassWriter, i: number, ink: number, rows = 1) {
  const { box } = lay;
  const y0 = Math.max(lay.rowY(i), i >= lay.display.pinned ? box.y + lay.pinnedHeight : box.y);
  const y1 = Math.min(lay.rowY(i) + rows * lay.rowScale, box.y + box.height);
  if (y1 <= y0) return;
  const g = out.args;
  g[0] = box.x;
  g[1] = y0;
  g[2] = box.width;
  g[3] = y1 - y0;
  g[4] = ((ink & 0xffffff00) | 0x1a) >>> 0;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/** The hovered cell's x extents (`hoverCell` reads them: set before the call). */
const CellX = new Float64Array(2);

/**
 * The hovered cell's outline at display row `i`, between `CellX[0]` and `CellX[1]`: a squircle by the cell's size (every emphasis a
 * squircle), its radius plus the offset, square on a side the view or the band's edge cuts and below the minimum.
 */
function hoverCell(lay: RowsLay, out: PassWriter, i: number, ink: number) {
  const { box, config } = lay;
  // the cell's x through an array: a fractional number passed to a call that is not inlined is boxed
  const a = CellX[0];
  const b = CellX[1];
  const y = lay.rowY(i);
  const rh = lay.rowScale;
  if (y + rh <= box.y || y >= box.y + box.height) return;
  const g = out.args;
  const d = config.outline.hover / 2 + config.emphasis.hoverOffset;
  const w = b - a;
  const r =
    w < config.squircleMin || rh < config.squircleMin
      ? 0
      : Math.min(rh * config.squircleRatio, config.squircleMax, w / 2, rh / 2);
  // cut: past the data box's sides (the view), above or below the band (its rows scroll inside it)
  let round = (a >= box.x ? 1 : 0) | (b <= box.x + box.width ? 2 : 0);
  if (y < box.y) round |= 4;
  if (y + rh > box.y + box.height) round |= 8;
  g[0] = a - d;
  g[1] = y - d;
  g[2] = w + 2 * d;
  g[3] = rh + 2 * d;
  g[4] = NoColor;
  g[5] = ink;
  g[6] = config.outline.hover;
  g[7] = r > 0 && (round & 3) !== 0 ? r + d : 0;
  g[8] = round;
  out.rectArgs();
}

/**
 * A hovered clade (_spec/README.md "Render", The rows layer, Hover): its branches in the hover ink, one band over each
 * run of its rows in the order (a folded clade's summary row alone), a pinned row's band too.
 */
function hoverClade(lay: RowsLay, out: PassWriter, key: string, ink: number) {
  const cl = lay.clades;
  const c = cl?.byName.get(key);
  if (!cl || c === undefined) return;
  const g = out.args;
  const { box, config } = lay;
  const P = lay.display.pinned;
  let start = -1;
  let end = -1;
  for (let l = cl.lo[c]; l < cl.hi[c]; l++) {
    const r = cl.leafRow[l];
    const p = r >= 0 ? lay.display.placeRow(r) : -1;
    if (p < 0) continue;
    // a pinned row in the band is lit there (a run is the scrolled rows and its place in the order)
    if (P) {
      const d = lay.display.displayRow(r);
      if (d >= 0 && d < P) hoverBand(lay, out, d, ink);
    }
    if (p < P) continue;
    if (start < 0) start = end = p;
    else if (p === end)
      continue; // a folded clade's rows share its summary row
    else if (p === end + 1) end = p;
    else {
      hoverBand(lay, out, start, ink, end - start + 1);
      start = end = p;
    }
  }
  if (start >= 0) hoverBand(lay, out, start, ink, end - start + 1);
  g[0] = config.inset;
  g[1] = box.y;
  g[2] = Math.max(0, box.x - config.inset);
  g[3] = box.height;
  out.clipArgs();
  const t = lay.tree;
  const owner = lay.treeClade;
  branchWidth(lay, config.outline.hover);
  for (let k = 0, s = 0; k + 3 < t.length; k += 4, s++) if (within(cl, owner[s], c)) branchRect(lay, out, k, ink);
  if (lay.wedges.length) drawWedges(lay, out, ink, c, true);
  out.unclip();
}

/**
 * The hover: its row's band tinted lightly (a pinned row at both its places: the pinned band and the order), its
 * column within the rows, the cell outlined.
 */
export function drawHover(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const h = emphasis?.hover;
  if (!h || h.layer !== lay.key) return;
  const { box, tokens } = lay;
  const g = out.args;
  const ink = tokens.colors['seq-hover'];
  if ((h.target === 'tree-node' || h.target === 'summary-row') && h.element?.owner === lay.owner)
    hoverClade(lay, out, String(h.element.key), ink);
  const row = h.yPosition ?? (h.system === lay.rowsSystem ? h.position : undefined);
  const i = row === undefined ? -1 : lay.display.displayRow(row);
  const j = row === undefined ? -1 : lay.display.placeRow(row);
  if (i >= 0) hoverBand(lay, out, i, ink);
  if (j >= 0 && j !== i) hoverBand(lay, out, j, ink);
  if (h.position === undefined || h.system !== lay.system || !lay.axis.cell(h.position)) return;
  const dpr = lay.axis.dpr;
  const a = Math.round(lay.axis.sx0 * dpr) / dpr;
  const b = Math.max(Math.round(lay.axis.sx1 * dpr) / dpr, a + 1 / dpr);
  g[0] = a;
  g[1] = box.y;
  g[2] = b - a;
  g[3] = box.height;
  g[4] = ((ink & 0xffffff00) | 0x14) >>> 0;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
  CellX[0] = a;
  CellX[1] = b;
  if (i >= 0) hoverCell(lay, out, i, ink);
  if (j >= 0 && j !== i) hoverCell(lay, out, j, ink);
}

/** The selection's tint and outline colors and the x span (`drawSelection` sets them): read by `paintRows`. */
const SelPaint = new Float64Array(4);

/**
 * The selection along y: a rows part (data rows) as row bands over the data box, tinted with an outline; a rect of
 * rows × columns as its cells' boxes. A summary row draws the partial mark: the tint over the share of its rows
 * selected, the outline when all are. O(rows in view) whatever the selection's size.
 */
export function drawSelection(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const sel = emphasis?.selection;
  if (!sel) return;
  const { box, tokens, config } = lay;
  const accent = tokens.colors['msui-accent'];
  // the colors through a typed array: a number past a small integer passed to a call is boxed
  SelPaint[0] = ((accent & 0xffffff00) | Math.round(255 * config.emphasis.selectionTint)) >>> 0;
  SelPaint[1] = ((accent & 0xffffff00) | Math.round(255 * config.emphasis.selectionOutline)) >>> 0;
  for (let p = 0; p < sel.parts.length; p++) {
    const part = sel.parts[p];
    if (part.system !== lay.rowsSystem) continue;
    SelPaint[2] = box.x;
    SelPaint[3] = box.x + box.width;
    paintRows(lay, out, part.ranges);
  }
  for (let r = 0; r < sel.rects.length; r++) {
    const rect = sel.rects[r];
    const rows = rect.a.system === lay.rowsSystem ? rect.a : rect.b.system === lay.rowsSystem ? rect.b : undefined;
    const cols = rows === rect.a ? rect.b : rect.a;
    if (!rows || cols.system !== lay.system) continue;
    for (let c = 0; c + 1 < cols.ranges.length; c += 2) {
      const xa = lay.axis.x(cols.ranges[c]);
      const xb = lay.axis.x(cols.ranges[c + 1]);
      if (Number.isNaN(xa) || Number.isNaN(xb)) continue;
      SelPaint[2] = Math.min(xa, xb);
      SelPaint[3] = Math.max(xa, xb);
      paintRows(lay, out, rows.ranges);
    }
  }
}

/** Whether a data row is in sorted half-open `ranges` (a binary search). */
function inRanges(ranges: ArrayLike<number>, row: number): boolean {
  let lo = 0;
  let hi = ranges.length >> 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (ranges[2 * mid + 1] <= row) lo = mid + 1;
    else hi = mid;
  }
  return lo < ranges.length >> 1 && ranges[2 * lo] <= row;
}

/** The display rows in view (the pinned band, then the scrolled rows) whose data rows `ranges` holds, each once. */
function paintRows(lay: RowsLay, out: PassWriter, ranges: ArrayLike<number>) {
  for (let i = 0; i < lay.display.pinned; i++) paintRow(lay, out, i, ranges);
  for (let i = lay.first; i < lay.last; i++) paintRow(lay, out, i, ranges);
}

function paintRow(lay: RowsLay, out: PassWriter, i: number, ranges: ArrayLike<number>) {
  const e = lay.display.rowAt(i);
  let share = 0;
  if (e >= 0) share = inRanges(ranges, e) ? 1 : 0;
  else {
    const m = lay.display.members.get(i);
    if (!m) return;
    let k = 0;
    for (let q = 0; q < m.length; q++) if (inRanges(ranges, m[q])) k++;
    share = k / m.length;
  }
  if (share <= 0) return;
  const { box } = lay;
  const y = lay.rowY(i);
  if (y + lay.rowScale <= box.y || y >= box.y + box.height) return;
  const g = out.args;
  g[0] = SelPaint[2];
  g[1] = y;
  g[2] = (SelPaint[3] - SelPaint[2]) * share;
  g[3] = lay.rowScale;
  g[4] = SelPaint[0];
  g[5] = share === 1 ? SelPaint[1] : NoColor;
  g[6] = share === 1 ? lay.config.outline.selection : 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}

/**
 * The keyboard cursor on a cell of two axes (_spec/README.md "Render", The rows layer), on its own section: its cell, or
 * its row's header (the label column, `labelX` to `numberX`), outlined in the neutral ink, a squircle by the cell (every
 * emphasis a squircle) inside it by half its width, square on a side the view or the band cuts and below the minimum. A
 * folded row's at its summary row, a hidden one's at the nearest row shown (`cursorRow`); a pinned row's in the band and
 * in its place too (by data row, as the selection).
 */
export function drawCursor(lay: RowsLay, out: PassWriter, emphasis: Emphasis | undefined) {
  const c = emphasis?.cursor;
  if (!c || c.ySystem !== lay.rowsSystem || c.yPosition === undefined) return;
  // the layer's address is its section's and more: `startsWith` and a code read allocate nothing
  const s = c.section;
  if (s !== undefined && !(lay.key.startsWith(s) && lay.key.charCodeAt(s.length) === 47)) return;
  const i = lay.cursorRow(c.yPosition);
  if (i < 0) return;
  const header = c.part === 'header';
  if (header) {
    CellX[0] = lay.labelX;
    CellX[1] = lay.numberX;
  } else {
    if (c.system !== lay.system || !lay.axis.cell(c.position)) return;
    CellX[0] = Math.min(lay.axis.sx0, lay.axis.sx1);
    CellX[1] = Math.max(lay.axis.sx0, lay.axis.sx1);
  }
  cursorCell(lay, out, i, header);
  const own = lay.display.placeRow(c.yPosition);
  if (own >= 0 && own !== i) cursorCell(lay, out, own, header);
}

/** The cursor's outline at display row `i` between `CellX[0]` and `CellX[1]` (`drawCursor`). */
function cursorCell(lay: RowsLay, out: PassWriter, i: number, header: boolean) {
  const { box, tokens, config } = lay;
  const a = CellX[0];
  const b = CellX[1];
  const y = lay.rowY(i);
  const rh = lay.rowScale;
  if (y + rh <= box.y || y >= box.y + box.height || !(b > a)) return;
  const w = config.outline.cursor;
  const cw = b - a;
  const r = squircleRadius(cw, rh, config.squircleRatio, config.squircleMax, config.squircleMin);
  // cut: past the data box's sides (the view; the header is whole), above or below the band (its rows scroll in it)
  let round = header ? Round.both : (a >= box.x ? Round.left : 0) | (b <= box.x + box.width ? Round.right : 0);
  if (y < box.y) round |= Round.cutTop;
  if (y + rh > box.y + box.height) round |= Round.cutBottom;
  const g = out.args;
  g[0] = a + w / 2;
  g[1] = y + w / 2;
  g[2] = Math.max(cw - w, 0.5);
  g[3] = Math.max(rh - w, 0.5);
  g[4] = NoColor;
  g[5] = tokens.colors[config.emphasis.cursorInk];
  g[6] = w;
  g[7] = r > 0 && (round & Round.both) !== 0 ? Math.max(r - w / 2, 0) : 0;
  g[8] = round;
  out.rectArgs();
}
