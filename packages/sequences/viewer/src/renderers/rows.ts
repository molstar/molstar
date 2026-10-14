/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `rows` renderer (seqspec:viewer#extension-points; seqspec:viewer#the-render-backend,
// seqspec:sources#levels-and-the-zoom-gate): an alignment's rows from tile images (renderers/rows-tiles.ts: held at the
// finest level, the region cache's reduced tiles, coarse under fine), the gaps between wide cells and the letters over
// them at the finest level (held codes or `tile(0, 0, …)`, in the ink that reads on their cell), each row's first and
// last residue numbers in the gutters, the whole alignment for the vertical overview, the tree in the header. Levels
// per block from `blockLevels` over `DrawView.lanes`, the same the stream manager requests. Batches are this module's,
// reset at a draw's start and handed over once at its end.

import { TILE, type Tile } from "@molstar/seq-core/format";
import { NO_LABEL } from "@molstar/seq-core/model";
import {
  type Axis,
  blockLevels,
  CHARS,
  cellInk,
  type DrawView,
  drawTree,
  drawTreeHover,
  kernelEntry,
  LETTER_LINE,
  letterFontPx,
  letterLines,
  letterRoom,
  NO_COLOR,
  numberText,
  overviewLevels,
  RectBatch,
  RowsHeader,
  type RowsState,
  rowAt,
  rowLabel,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
  ThemeLocation,
  textShown,
  wedgeLabel,
} from "../draw/api.ts";
import { beginTiles, paintLevels, paintTiles, RowTiles, rowGeo } from "./rows-tiles.ts";

const gaps = new RectBatch(256);
/** Insert columns' cells, veiled in the surface's tone. */
const muted = new RectBatch(64);
/** The veil's opacity (of 255). */
const MUTE = 140;
const letters = new TextBatch(256);
const numbers = new TextBatch(2);
const lv = new Int32Array(4);
/** The baselines that center a capital and a lower-case letter on their cell (`letterLines`). */
const baselines = new Float64Array(2);

/** A display row's lane offset (the layout's `RowModel`); -1: a collapsed clade's wedge. */
const laneAt = (r: NonNullable<Axis["rows"]>, d: number) =>
  r.node && r.node[d] >= 0 ? -1 : r.lane ? r.lane[d] - r.laneStart : d;
/** Whether block b shows this layer's alignment and is drawn here. */
const shows = (v: DrawView, s: RowsState, b: number) =>
  (v.axis as Axis).alignment[b] === s.alignment && (v.onlyBlock < 0 || v.onlyBlock === b) && v.visHi[b] > v.visLo[b];

/** A draw's start: the cache's images checked, the batches cleared. */
function begin(v: DrawView, c: RowTiles): void {
  beginTiles(v, c);
  gaps.reset();
  muted.reset();
  letters.reset();
  ThemeLocation.reset(v.loc);
}
/** A draw's end: the insert columns' veil and the gaps over the cells, the letters over them. */
function end(v: DrawView): void {
  if (muted.n > 0) v.backend.rects(muted);
  if (gaps.n > 0) v.backend.rects(gaps);
  if (letters.n > 0) v.backend.text(letters);
}

/**
 * Lanes [k, k + n) of every block that shows them: the tiles the block's levels ask for, at `rowGeo`; `open` as
 * `paintTiles`' (the lanes above, below, outside the clip).
 */
function cells(v: DrawView, s: RowsState, c: RowTiles, k: number, n: number, open = 0): void {
  const budget = v.stream ? v.stream.budget : 0;
  for (let b = 0; b < (v.axis as Axis).count; b++) {
    if (!shows(v, s, b)) continue;
    blockLevels(v, b, budget, v.lanes, lv, v.stream?.served);
    paintLevels(v, s, c, b, lv, k, n, v.visLo[b], v.visHi[b], open);
  }
}

/** One row (lane offset `k`, at `rowGeo`): its cells, the gaps between wide cells, its letters. */
function row(v: DrawView, s: RowsState, c: RowTiles, k: number): void {
  const ax = v.axis as Axis;
  const y0 = Math.round(rowGeo[0] * v.dpr) | 0;
  const y1 = Math.round((rowGeo[0] + rowGeo[1]) * v.dpr) | 0;
  const slot = s.slot ? s.slot[k] : -1;
  cells(v, s, c, k, 1);
  // The row's text as the letters decide it (`textShown`); the letter's size from the row's height.
  letterRoom[0] = rowGeo[1];
  letterRoom[1] = 0;
  for (let b = 0; b < ax.count; b++) {
    if (!shows(v, s, b)) continue;
    const gap = v.cell * v.span[b] > 6 ? 1 : 0;
    const lo = v.visLo[b];
    const hi = v.visHi[b];
    if (gap > 0) {
      gaps.begin(SkinToken.surface);
      for (let q = lo; q < hi; q++) gaps.push(v.xe(b, q) - gap, y0, gap, y1 - y0);
    }
    // Insert columns (seqspec:model#alignments-and-trees): their cells veiled, their letters lower case.
    const ins = ax.alignment[b] >= 0 ? (ax.set.alignments[ax.alignment[b]]?.insertColumn ?? null) : null;
    if (ins !== null) {
      muted.begin(SkinToken.surface, MUTE);
      for (let q = lo; q < hi; q++)
        if (ins[ax.indexOf(b, q as never)] === 1) muted.push(v.xd(b, q), y0, v.xe(b, q) - v.xd(b, q), y1 - y0);
    }
    if (!textShown(v, b)) continue;
    const font = v.text.font("mono", letterFontPx(v, b));
    if (letters.n > 0 && letters.font !== font) {
      v.backend.text(letters);
      letters.reset();
    }
    letters.style(font, TextAlign.center, TextBaseline.alphabetic);
    // On the cell `paintTiles` draws: the row inset by `rowGeo[2]` above and below.
    const inset = rowGeo[2];
    letterLines(
      v,
      font,
      Math.round((rowGeo[0] + inset) * v.dpr) | 0,
      Math.round((rowGeo[0] + rowGeo[1] - inset) * v.dpr) | 0,
      baselines,
    );
    let tile: Tile | null = null;
    let tx = -1;
    const ar = rowAt(s, k);
    for (let q = lo; q < hi; q++) {
      const i = ax.indexOf(b, q as never);
      let code = 0;
      if (s.codes) code = slot >= 0 && i < s.columns ? s.codes[slot * s.columns + i] : 0;
      else if (v.stream) {
        if (i >> 8 !== tx) {
          tx = i >> 8;
          tile = v.stream.tile(0, 0, k >> 8, tx);
        }
        code = tile ? tile.cls[(k & 255) * TILE + (i & 255)] : 0;
      }
      if (code === 0) continue;
      if (s.map !== null) code = s.map[code & 127];
      if (ins !== null && ins[i] === 1 && code >= 65 && code <= 90) code += 32;
      let fill = NO_COLOR;
      if (v.color) {
        const loc = ThemeLocation.at(v.loc, ax, v.layer, b, i);
        loc.row = ar;
        loc.letter = code;
        fill = v.color.value(loc);
      }
      letters.begin(code === 45 || code === 46 ? SkinToken.inkQuiet : fill >= 0 ? cellInk(fill) : SkinToken.ink);
      const t = letters.add(CHARS[code & 127]);
      letters.x[t] = (v.xd(b, q) + v.xe(b, q) - gap) / 2;
      letters.y[t] = baselines[LETTER_LINE[code & 127]];
    }
  }
}

/** One row (`v.row`, a lane offset) in the entry's box. */
function draw(v: DrawView, s: RowsState): void {
  const c = v.cache as RowTiles | null;
  if (!c || v.row < 0 || !v.axis || v.lanes[1] <= v.lanes[0]) return;
  begin(v, c);
  rowGeo[0] = v.top;
  rowGeo[1] = v.height;
  rowGeo[2] = v.height >= 8 ? 1 : 0;
  row(v, s, c, v.row);
  end(v);
}

/**
 * Display rows [first, last] in one call (`v.top`: row `first`'s): row by row at a readable height, else as runs of
 * contiguous lanes, one image draw per tile.
 */
function rows(v: DrawView, s: RowsState, _p: unknown, first: number, last: number): void {
  const r = v.axis?.rows;
  const c = v.cache as RowTiles | null;
  if (!r || !c || v.lanes[1] <= v.lanes[0]) return;
  begin(v, c);
  const dense = v.height < v.props.denseRowHeight;
  for (let d = first; d <= last; ) {
    const k = laneAt(r, d);
    if (k < 0) {
      d++;
      continue;
    }
    rowGeo[0] = v.top + (d - first) * v.height;
    rowGeo[1] = v.height;
    rowGeo[2] = !dense && v.height >= 8 ? 1 : 0;
    if (!dense) {
      row(v, s, c, k);
      d++;
      continue;
    }
    let e = d + 1;
    while (e <= last && laneAt(r, e) === k + (e - d)) e++;
    // The call's first and last rows are the clip's (the composer's band): past them, nothing shows.
    cells(v, s, c, k, e - d, (d === first ? 1 : 0) | (e - 1 === last ? 2 : 0));
    d = e;
  }
  end(v);
}

/**
 * The whole alignment into the view's box (`v` a map, draw/miniature.ts' `MapView`): at the level the region cache
 * loads for the section's overview column (`overviewLevels`, whatever the box), the coarsest holding every row and
 * column in a tile or two per axis; a tile not yet arrived is the skin's `pending` shade, never the ground. Without an
 * overview column, at `MAP_WIDTH`'s (the manager loads them where a miniature names the track).
 */
function overview(v: DrawView, s: RowsState): boolean {
  const r = v.axis?.rows;
  const c = v.cache as RowTiles | null;
  if (!r || !c || r.count === 0 || s.columns === 0 || !shows(v, s, 0)) return false;
  begin(v, c);
  overviewLevels(v, r.laneCount, s.columns, v.stream ? v.stream.budget : 0, lv, v.stream?.served);
  rowGeo[1] = v.height / r.count;
  rowGeo[2] = 0;
  for (let d = 0; d < r.count; ) {
    const k = laneAt(r, d);
    let e = d + 1;
    while (k >= 0 && e < r.count && laneAt(r, e) === k + (e - d)) e++;
    rowGeo[0] = v.top + d * rowGeo[1];
    if (k >= 0) paintTiles(v, s, c, 0, lv[2], lv[3], k, e - d, v.visLo[0], v.visHi[0], SkinToken.pending);
    d = e;
  }
  return true;
}

/**
 * A row's first and last residue numbers in view, in the left and right gutters (rows held in memory). An alignment
 * axis is one alignment's columns, one block (`svs.section.axis`): block 0.
 */
function gutter(v: DrawView, s: RowsState): void {
  const a = s.set.alignments[s.alignment];
  if (!a || !v.axis || a.lazy || v.dataRow < 0 || v.height < 10 || v.gutterLeft < 20 || !shows(v, s, 0)) return;
  const slot = a.rows.slotOf(v.dataRow as never);
  if (slot < 0) return;
  const first = a.rows.rowFirst[slot] === NO_LABEL ? 1 : a.rows.rowFirst[slot];
  // Digits centered on their height in the row, as the row's letters (`letterLines`; a `middle` baseline drew them high).
  const font = v.text.font("ui", 10);
  letterLines(v, font, Math.round(v.top * v.dpr) | 0, Math.round((v.top + v.height) * v.dpr) | 0, baselines);
  const y = baselines[0];
  for (let side = 0; side < 2; side++) {
    let p = -1;
    for (
      let q = side === 0 ? v.visLo[0] : v.visHi[0] - 1;
      p < 0 && q >= v.visLo[0] && q < v.visHi[0];
      q += 1 - 2 * side
    )
      p = a.rows.colToPos(v.dataRow as never, v.axis.indexOf(0, q as never) as never);
    if (p < 0) continue;
    numbers.reset();
    numbers.style(font, side === 0 ? TextAlign.right : TextAlign.left, TextBaseline.alphabetic);
    numbers.begin(SkinToken.inkQuiet);
    const t = numbers.add(numberText(first + p));
    numbers.x[t] = Math.round((v.left - 4 + side * (v.right - v.left + 8)) * v.dpr);
    numbers.y[t] = y;
    v.backend.text(numbers);
  }
}

export const RowsRenderer = kernelEntry<RowsState>("rows", (_s, _p, props) => props.rowHeight, {
  minHeight: 1,
  miniature: false,
  draw,
  tiny: draw,
  rows,
  gutter,
  overview,
  header: {
    ...RowsHeader,
    draw: (v, s, _p, first, last) => drawTree(v, s, first, last),
    hover: (v, s, _p, node, first, last) => drawTreeHover(v, s, node, first, last),
  },
  rowModel: { label: rowLabel, wedge: wedgeLabel },
  cache: { create: () => new RowTiles(), dispose: (c) => (c as RowTiles).dispose() },
});
