/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { fontSize } from './config.js';
import { Anchor, NoColor, type PassWriter, Round } from './draw-lists.js';
import { RenderIssues } from './issues.js';
import type { LayoutContext } from './layer-type.js';
import type { RowsLay, RowsPrep } from './rows.js';
import { columnsRoom, formatColumnValue, numericValues, type RowColumn, RowMeasureFields } from './rows-model.js';

// The rows layer's row columns (_spec/README.md "Render", The rows layer, Row columns): each column's values resolved
// once a data (its field's origin, numbers where every value reads as one, its range and categories), a cell's text,
// its cut and its widths kept by row and room, the cells in view, the titles in the layer's own title row and the
// hidden-columns chip; drawn and hit from those (a layout once warm, a draw and a hover allocate nothing).

/** The tip's surface past its text on each side, px. */
const TipPad = 3;

const AnchorStart: number = Anchor.start;
const AnchorEnd: number = Anchor.end;

/** Where a column's field comes from: what formats it (a measure's `identity` is a percent, a joined one is not). */
export type FieldOrigin = 'label' | 'first' | 'own' | 'measure' | 'joined' | 'none';

/** The order marks a title carries when the rows are ordered by its field. */
const OrderMark = { ascending: ' ↑', descending: ' ↓' } as const;

/** The categories' colors a chip takes, by its category's index (the `categories` list). */
const CategoryTokens = [
  'seq-cat-0',
  'seq-cat-1',
  'seq-cat-2',
  'seq-cat-3',
  'seq-cat-4',
  'seq-cat-5',
  'seq-cat-6',
  'seq-cat-7',
] as const;

const NoCuts: (string | undefined)[] = [];

/** A field's origin on a layer: `label`, `first`, the rows' own, a measure, a joined field, else none. */
export function fieldOrigin(lay: RowsLay, field: string): FieldOrigin {
  if (field === 'label') return 'label';
  if (field === 'first') return 'first';
  if (lay.held?.fields?.has(field)) return 'own';
  if (lay.held && RowMeasureFields.includes(field)) return 'measure';
  if (lay.joined?.has(field)) return 'joined';
  return 'none';
}

/** `text` cut to `max` px with an ellipsis, as `font` measures it (whole when it fits). */
export function elide(text: LayoutContext['text'], font: string, t: string, max: number): string {
  if (text.measure(font, t) <= max) return t;
  let lo = 0;
  let hi = t.length;
  // the longest prefix that fits with its ellipsis
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (text.measure(font, `${t.slice(0, mid)}…`) <= max) lo = mid;
    else hi = mid - 1;
  }
  return lo ? `${t.slice(0, lo)}…` : '';
}

/** A chip's swatch, px: a square of the row's middle. */
function chipSwatch(rs: number): number {
  return Math.min(10, rs * 0.5);
}

/** What a rows layer keeps of its row columns, frame to frame. */
export class ColumnsLay {
  /** The columns as read, their widths after the view's yielding (0: hidden), their x, the count hidden. */
  columns: readonly RowColumn[] = [];
  width: number[] = [];
  x: number[] = [];
  hidden = 0;
  /** Where the right columns end (their clip). */
  rightEnd = 0;
  /** The layer's own title row above its rows: its height, and the titles' baseline in it. */
  titleRow = 0;
  titleY = 0;
  /**
   * Each column's field resolved over the held rows (kept by `of`): its origin, a value a data row (`null` missing;
   * numbers when every value is one or reads as one), whether it is numeric, its range (a bar's domain, from 0), its
   * categories in code-point order (a chip's color) and whether any row has it.
   */
  origin: FieldOrigin[] = [];
  values: (string | number | null)[][] = [];
  numeric: boolean[] = [];
  lo: number[] = [];
  hi: number[] = [];
  cats: Map<string, number>[] = [];
  any: boolean[] = [];
  ofHeld: unknown;
  ofFields = -1;
  ofPrep: unknown;
  /**
   * A cell's text by column and data row (made when first shown, kept with the values), what a hover names, its cut
   * at the column's room (`cutRoom`: a new room cuts again) and the two widths (NaN until measured).
   */
  full: (string | undefined)[][] = [];
  name: (string | undefined)[][] = [];
  cut: (string | undefined)[][] = [];
  cutRoom: number[] = [];
  cutW: Float64Array[] = [];
  fullW: Float64Array[] = [];
  /**
   * A cut cell's text as the hover reads it (the `tip` pass): whole, or, when the view's edge is nearer than its end,
   * cut there with an ellipsis; its width, and the room (to the view's edge from the text's x) they were cut for.
   */
  tip: (string | undefined)[][] = [];
  tipW: Float64Array[] = [];
  tipRoom: number[] = [];
  /** The cells in view: their column, display row, data row, baseline and row top, bar share, chip color; `n`. */
  n = 0;
  col = new Int32Array(64);
  display = new Int32Array(64);
  row = new Int32Array(64);
  y = new Float64Array(64);
  top = new Float64Array(64);
  share = new Float64Array(64);
  color = new Uint32Array(64);
  /** The titles: text, left x, drawn width, column; `titleN`; kept by what they were made for. */
  titleN = 0;
  titleText: string[] = [];
  titleX: number[] = [];
  titleW: number[] = [];
  titleCol: number[] = [];
  titlesFor: unknown[] = [undefined, undefined, Number.NaN, Number.NaN];
  /** The hidden-columns chip in the title row: its text, x and width ('' none). */
  chip = '';
  chipX = 0;
  chipW = 0;
  /** A field no row has (`seq.render.rows.column`), made once a data, said every frame it stands. */
  missing: Issue | undefined;
  missingOf = -1;
  /** The font's numbers and size, measured once a font. */
  font = '';
  numberW = 0;
  size = 0;

  /** The text a cell of column `k` shows at data row `r` whole, made and kept on first use. */
  fullText(k: number, r: number): string {
    const c = this.full[k];
    let t = c[r];
    if (t === undefined) {
      const v = this.values[k][r];
      t = v === null ? '' : formatColumnValue(this.origin[k] === 'measure' ? this.columns[k].field : '', v);
      c[r] = t;
    }
    return t;
  }

  /** What a hover on the cell names: `title: value`, kept. */
  nameOf(k: number, r: number): string {
    const c = this.name[k];
    let t = c[r];
    if (t === undefined) c[r] = t = `${this.columns[k].title}: ${this.fullText(k, r)}`;
    return t;
  }

  /** The cell at display row `i` under `x`: its index in the cells, -1 for none. */
  cellAt(i: number, x: number): number {
    for (let q = 0; q < this.n; q++) {
      if (this.display[q] !== i) continue;
      const k = this.col[q];
      const x0 = this.x[k];
      if (x >= x0 && x < x0 + this.width[k]) return q;
    }
    return -1;
  }

  /** The cut text of cell `q`. */
  cellText(q: number): string {
    return this.cut[this.col[q]][this.row[q]] ?? '';
  }

  /** The whole text of cell `q`. */
  cellFull(q: number): string {
    return this.full[this.col[q]][this.row[q]] ?? '';
  }
}

/** Sets the columns as written and their granted widths, and the layer's title row (each frame; nothing allocated). */
export function setColumns(cols: ColumnsLay, prep: RowsPrep, ctx: LayoutContext) {
  const columns = prep.columns;
  if (cols.columns !== columns) {
    cols.columns = columns;
    cols.width = new Array(columns.length).fill(0);
    cols.x = new Array(columns.length).fill(0);
    cols.ofPrep = undefined;
  }
  const granted = ctx.rowColumns?.widths;
  cols.hidden = 0;
  for (let k = 0; k < columns.length; k++) {
    cols.width[k] = granted && granted.length === columns.length ? granted[k] : columns[k].width;
    if (cols.width[k] <= 0) cols.hidden++;
  }
  cols.titleRow = columns.length ? ctx.config.rows.columns.titlePx : 0;
  cols.n = 0;
}

/**
 * Each column's field resolved over the held rows: kept while the held data, the joined fields and the columns are
 * the same (one pass a column when one changes).
 */
function resolveColumns(lay: RowsLay, cols: ColumnsLay) {
  const held = lay.held;
  if (cols.ofHeld === held && cols.ofFields === lay.fieldsVersion && cols.ofPrep === cols.columns) return;
  cols.ofHeld = held;
  cols.ofFields = lay.fieldsVersion;
  cols.ofPrep = cols.columns;
  cols.missing = undefined;
  cols.missingOf = -1;
  const n = held?.rows ?? 0;
  const columns = cols.columns;
  for (let k = 0; k < columns.length; k++) {
    const field = columns[k].field;
    cols.origin[k] = fieldOrigin(lay, field);
    const raw: (string | number | null)[] = new Array(n);
    for (let r = 0; r < n; r++) raw[r] = lay.fieldValue(r, field) ?? null;
    // a column whose values all read as numbers (a `#=GS` TAXID, an untyped table's field) is numeric for its purposes
    const nums = numericValues(raw);
    const values = nums ?? raw;
    let lo = 0;
    let hi = Number.NEGATIVE_INFINITY;
    let any = false;
    const seen = new Set<string>();
    for (let r = 0; r < n; r++) {
      const v = values[r];
      if (v === null) continue;
      any = true;
      if (typeof v === 'number') {
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      } else seen.add(v);
    }
    cols.values[k] = values;
    cols.numeric[k] = nums !== undefined;
    cols.lo[k] = lo;
    cols.hi[k] = Number.isFinite(hi) ? hi : lo;
    cols.any[k] = any;
    const cats = new Map<string, number>();
    for (const c of [...seen].sort()) cats.set(c, cats.size);
    cols.cats[k] = cats;
    cols.full[k] = new Array(n);
    cols.name[k] = new Array(n);
    cols.cut[k] = NoCuts;
    cols.cutRoom[k] = Number.NaN;
    cols.cutW[k] = new Float64Array(n).fill(Number.NaN);
    cols.fullW[k] = new Float64Array(n).fill(Number.NaN);
    cols.tip[k] = NoCuts;
    cols.tipRoom[k] = Number.NaN;
    cols.tipW[k] = new Float64Array(n).fill(Number.NaN);
  }
  for (const a of [cols.origin, cols.values, cols.numeric, cols.lo, cols.hi, cols.cats, cols.any])
    a.length = columns.length;
  cols.full.length = cols.name.length = cols.cut.length = columns.length;
  cols.cutRoom.length = cols.cutW.length = cols.fullW.length = columns.length;
  cols.tip.length = cols.tipW.length = cols.tipRoom.length = columns.length;
}

function grow(cols: ColumnsLay, need: number) {
  if (need <= cols.col.length) return;
  const m = Math.max(need, 2 * cols.col.length);
  const i32 = (a: Int32Array) => {
    const b = new Int32Array(m);
    b.set(a);
    return b;
  };
  const f64 = (a: Float64Array) => {
    const b = new Float64Array(m);
    b.set(a);
    return b;
  };
  cols.col = i32(cols.col);
  cols.display = i32(cols.display);
  cols.row = i32(cols.row);
  cols.y = f64(cols.y);
  cols.top = f64(cols.top);
  cols.share = f64(cols.share);
  const c = new Uint32Array(m);
  c.set(cols.color);
  cols.color = c;
}

/** The titles in the title row and the hidden-columns chip, made again only when what they read changed. */
function layoutTitles(lay: RowsLay, cols: ColumnsLay, prep: RowsPrep, ctx: LayoutContext, end: number) {
  const f = cols.titlesFor;
  const granted = ctx.rowColumns?.widths ?? prep.columns;
  if (f[0] === prep && f[1] === granted && f[2] === end && f[3] === cols.x[0]) return;
  f[0] = prep;
  f[1] = granted;
  f[2] = end;
  f[3] = cols.x[0];
  const { text, config } = ctx;
  const font = config.itemFont;
  const gap = config.rows.columns.gapPx;
  const columns = cols.columns;
  // the hidden-columns chip ends where the right columns do (before a rows rail): the right titles cut for it
  let limit = Number.POSITIVE_INFINITY;
  cols.chip = '';
  if (cols.hidden) {
    cols.chip = cols.hidden === 1 ? '1 column hidden' : `${cols.hidden} columns hidden`;
    cols.chipW = text.measure(font, cols.chip);
    cols.chipX = end - cols.chipW;
    limit = cols.chipX - gap;
  }
  const o = prep.order;
  let t = 0;
  for (let k = 0; k < columns.length; k++) {
    const c = columns[k];
    const w = c.side === 'right' ? Math.min(cols.width[k], limit - cols.x[k]) : cols.width[k];
    if (w <= 0) continue;
    const mark = typeof o === 'object' && o.field === c.field ? OrderMark[o.direction] : '';
    const s = elide(text, font, `${c.title}${mark}`, w);
    if (!s) continue;
    const tw = text.measure(font, s);
    cols.titleText[t] = s;
    cols.titleW[t] = tw;
    cols.titleX[t] = c.as === 'number' ? cols.x[k] + w - tw : cols.x[k];
    cols.titleCol[t] = k;
    t++;
  }
  cols.titleN = t;
}

/**
 * The row columns laid out (_spec/README.md "Render", The rows layer, Row columns): the left ones in the header between
 * the label and the first numbers, the right ones after the last numbers (before a rows rail); the titles in the layer's
 * title row (a figure's first line only); a cell a row in view (a summary row's empty, a missing value empty). Once
 * warm (the columns resolved, the rows' texts cut) it allocates nothing.
 */
export function layoutRowColumns(lay: RowsLay, prep: RowsPrep, ctx: LayoutContext) {
  const cols = lay.cols;
  const { config, text } = ctx;
  const box = lay.box;
  const columns = cols.columns;
  const gap = config.rows.columns.gapPx;
  const font = config.itemFont;
  if (cols.font !== font) {
    cols.font = font;
    cols.numberW = text.measure(font, '00000');
    cols.size = fontSize(font);
  }
  const numberW = prep.numbers === 'gutters' ? cols.numberW : 0;
  // left: ending where the first numbers start
  let x = box.x - config.labelGap - numberW - 4 - columnsRoom(columns, cols.width, 'left', gap);
  for (let k = 0; k < columns.length; k++)
    if (columns[k].side === 'left' && cols.width[k] > 0) {
      cols.x[k] = x + gap;
      x += gap + cols.width[k];
    }
  // right: after the last numbers' gutter, where the view puts them (a figure's: after its lines)
  let r = ctx.rowColumns?.rightX ?? box.x + box.width + (numberW && !lay.figure ? numberW + config.labelGap : 0);
  for (let k = 0; k < columns.length; k++)
    if (columns[k].side === 'right' && cols.width[k] > 0) {
      cols.x[k] = r + gap;
      r += gap + cols.width[k];
    }
  cols.rightEnd = r;
  resolveColumns(lay, cols);
  if (lay.held && !lay.joinPending) {
    if (cols.missingOf !== lay.fieldsVersion) {
      cols.missingOf = lay.fieldsVersion;
      let names = '';
      for (let k = 0; k < columns.length; k++)
        if (!cols.any[k]) names += names ? `, ${columns[k].field}` : columns[k].field;
      cols.missing = names
        ? issueOf(RenderIssues, 'seq.render.rows.column', `No row has the field ${names}`, {
            address: ctx.key,
            hint: RenderIssues.codes['seq.render.rows.column'].hint,
          })
        : undefined;
    }
    if (cols.missing) ctx.report(cols.missing);
  }
  cols.titleY = box.y - 4;
  if (ctx.strip === undefined || ctx.strip === 0) layoutTitles(lay, cols, prep, ctx, r);
  else cols.titleN = 0;
  const rs = lay.rowScale;
  if (rs < config.rows.labelMinPx) return;
  for (let i = 0; i < lay.display.pinned; i++) addCells(lay, cols, ctx, i);
  for (let i = lay.first; i < lay.last; i++) addCells(lay, cols, ctx, i);
}

/**
 * A cut cell's text for the hover, made once a row and room: the whole text when it ends before the view's edge, else cut
 * with an ellipsis where the edge leaves room (a number's grows left, to its column's x). The room a column's text has to
 * the edge is its x past the swatch (`inset`) to the view's width, less the surface's padding.
 */
function tipOf(cols: ColumnsLay, k: number, e: number, ctx: LayoutContext, inset: number) {
  const right = cols.columns[k].as === 'number';
  const edge = ctx.viewWidth ?? Number.POSITIVE_INFINITY;
  const room = right ? cols.x[k] + cols.width[k] - TipPad : edge - (cols.x[k] + inset) - TipPad;
  // a new room (another view width) makes every row's again
  if (cols.tipRoom[k] !== room) {
    cols.tipRoom[k] = room;
    cols.tip[k] = new Array(cols.values[k].length);
    cols.tipW[k].fill(Number.NaN);
  }
  if (cols.tip[k][e] !== undefined) return;
  const s = elide(ctx.text, ctx.config.itemFont, cols.full[k][e]!, room);
  cols.tip[k][e] = s;
  cols.tipW[k][e] = ctx.text.measure(ctx.config.itemFont, s);
}

/** The cells of display row `i`: their text cut at the column's room and measured, once a row and room. */
function addCells(lay: RowsLay, cols: ColumnsLay, ctx: LayoutContext, i: number) {
  const e = lay.display.rowAt(i);
  if (e < 0 || !lay.held) return;
  const { text, config, tokens } = ctx;
  const font = config.itemFont;
  const rs = lay.rowScale;
  const columns = cols.columns;
  const top = lay.rowY(i);
  const y = top + rs / 2 + cols.size * 0.35;
  for (let k = 0; k < columns.length; k++) {
    const w = cols.width[k];
    if (w <= 0) continue;
    const v = cols.values[k][e];
    if (v === null || v === undefined) continue;
    const c = columns[k];
    let share = Number.NaN;
    let color = NoColor;
    let room = w;
    if (c.as === 'bar') {
      const lo = cols.lo[k];
      const hi = cols.hi[k];
      share = typeof v === 'number' ? (hi > lo ? (v - lo) / (hi - lo) : 1) : Number.NaN;
      room = 0;
    } else if (c.as === 'chip') {
      const at = cols.cats[k].get(typeof v === 'string' ? v : cols.fullText(k, e)) ?? 0;
      color = tokens.colors[CategoryTokens[at % CategoryTokens.length]];
      room = Math.max(0, w - chipSwatch(rs) - 4);
    }
    if (room > 0) {
      // the cut is kept at its room: a new room (a narrower view) cuts every row again
      if (cols.cutRoom[k] !== room) {
        cols.cutRoom[k] = room;
        cols.cut[k] = new Array(cols.values[k].length);
        cols.cutW[k].fill(Number.NaN);
      }
      if (cols.cut[k][e] === undefined) {
        const full = cols.fullText(k, e);
        const s = elide(text, font, full, room);
        cols.cut[k][e] = s;
        cols.cutW[k][e] = text.measure(font, s);
        cols.fullW[k][e] = s === full ? cols.cutW[k][e] : text.measure(font, full);
      }
      if (cols.cut[k][e] !== cols.full[k][e]) tipOf(cols, k, e, ctx, c.as === 'chip' ? chipSwatch(rs) + 4 : 0);
      cols.nameOf(k, e);
    } else cols.nameOf(k, e);
    const q = cols.n++;
    if (q >= cols.col.length) grow(cols, q + 1);
    cols.col[q] = k;
    cols.display[q] = i;
    cols.row[q] = e;
    cols.y[q] = y;
    cols.top[q] = top;
    cols.share[q] = share;
    cols.color[q] = color;
  }
}

// ---- drawing (nothing allocated: numbers through the writer's `args`) ----

/**
 * The cells of one side (0 left, 1 right): text from the column's x (a number's to its right end), a bar its share of
 * the width in the accent, a chip a swatch of its category's color before its text.
 */
export function drawColumnCells(lay: RowsLay, out: PassWriter, side: 0 | 1) {
  const { config, tokens } = lay;
  const cols = lay.cols;
  const g = out.args;
  const font = config.itemFont;
  const ink = tokens.colors['msui-text'];
  const rs = lay.rowScale;
  for (let q = 0; q < cols.n; q++) {
    const k = cols.col[q];
    const c = cols.columns[k];
    if ((c.side === 'right' ? 1 : 0) !== side) continue;
    const x = cols.x[k];
    const w = cols.width[k];
    if (c.as === 'bar') {
      const share = cols.share[q];
      if (Number.isNaN(share)) continue;
      g[0] = x;
      g[1] = cols.top[q] + rs * 0.25;
      g[2] = Math.max(1, w * share);
      g[3] = rs * 0.5;
      g[4] = tokens.colors['msui-accent'];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
      continue;
    }
    let tx = x;
    if (c.as === 'chip') {
      const sw = chipSwatch(rs);
      g[0] = x;
      g[1] = cols.top[q] + (rs - sw) / 2;
      g[2] = sw;
      g[3] = sw;
      g[4] = cols.color[q];
      g[5] = NoColor;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      out.rectArgs();
      tx = x + sw + 4;
    }
    const r = cols.row[q];
    const t = cols.cut[k][r];
    if (!t) continue;
    const number = c.as === 'number';
    g[0] = number ? x + w : tx;
    g[1] = cols.y[q];
    g[2] = cols.cutW[k][r];
    g[3] = cols.cutW[k][r];
    g[4] = ink;
    out.textArgs(t, number ? AnchorEnd : AnchorStart, font);
  }
}

/** The titles in the title row (muted), and the hidden-columns chip (the accent). */
export function drawTitles(lay: RowsLay, out: PassWriter) {
  const { config, tokens } = lay;
  const cols = lay.cols;
  const g = out.args;
  for (let t = 0; t < cols.titleN; t++) {
    g[0] = cols.titleX[t];
    g[1] = cols.titleY;
    g[2] = cols.titleW[t];
    g[3] = cols.titleW[t];
    g[4] = tokens.colors['msui-text-muted'];
    out.textArgs(cols.titleText[t], AnchorStart, config.itemFont);
  }
  if (cols.chip) {
    g[0] = cols.chipX;
    g[1] = cols.titleY;
    g[2] = cols.chipW;
    g[3] = cols.chipW;
    g[4] = tokens.colors['msui-accent'];
    out.textArgs(cols.chip, AnchorStart, config.itemFont);
  }
}

/**
 * A hovered display row's cut cells read whole, in the `tip` pass over the labels: each on an opaque surface (the page's
 * background, then the track's well color, outlined in the hover ink, a squircle by the shape rule) from the cell's x,
 * to the right up to the view's edge, where a longer text ends in an ellipsis. Nothing under it shows through.
 */
export function hoverColumnCells(lay: RowsLay, out: PassWriter, i: number) {
  const { tokens, config } = lay;
  const cols = lay.cols;
  const g = out.args;
  for (let q = 0; q < cols.n; q++) {
    if (cols.display[q] !== i) continue;
    const k = cols.col[q];
    const r = cols.row[q];
    const cut = cols.cut[k][r];
    const full = cols.full[k][r];
    if (!cut || full === undefined || cut === full) continue;
    const text = cols.tip[k][r];
    if (!text) continue;
    const c = cols.columns[k];
    const w = cols.tipW[k][r];
    const x =
      c.as === 'number'
        ? cols.x[k] + cols.width[k] - w
        : cols.x[k] + (c.as === 'chip' ? chipSwatch(lay.rowScale) + 4 : 0);
    const h = lay.rowScale;
    const sw = w + 2 * TipPad;
    const radius =
      sw < config.squircleMin || h < config.squircleMin
        ? 0
        : Math.min(h * config.squircleRatio, config.squircleMax, sw / 2, h / 2);
    for (let layer = 0; layer < 2; layer++) {
      g[0] = x - TipPad;
      g[1] = cols.top[q];
      g[2] = sw;
      g[3] = h;
      g[4] = layer === 0 ? tokens.colors['msui-bg'] : tokens.colors['seq-track'];
      g[5] = layer === 0 ? NoColor : tokens.colors['seq-hover'];
      g[6] = layer === 0 ? 0 : config.outline.hover;
      g[7] = radius;
      g[8] = radius > 0 ? Round.both : 0;
      out.rectArgs();
    }
    g[0] = x;
    g[1] = cols.y[q];
    g[2] = w;
    g[3] = w;
    g[4] = tokens.colors['msui-text'];
    out.textArgs(text, AnchorStart, config.itemFont);
  }
}

/** The title under a point: its index, -1 for none (its drawn text, in the title row; the chip is none). */
export function titleAt(lay: RowsLay, x: number, y: number): number {
  const cols = lay.cols;
  if (!cols.titleN) return -1;
  const top = lay.box.y - cols.titleRow;
  if (y < top || y >= lay.box.y) return -1;
  for (let t = 0; t < cols.titleN; t++) if (x >= cols.titleX[t] && x < cols.titleX[t] + cols.titleW[t]) return t;
  return -1;
}
