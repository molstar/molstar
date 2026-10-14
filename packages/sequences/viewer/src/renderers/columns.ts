/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The built-in data columns (seqspec:viewer#extension-points): `number` (right-aligned), `fraction` (a bar and a
// number), `text` (one line with an ellipsis), `category` (a swatch by the categorical palette and a label). `prepare`
// formats the strings once (seqspec:viewer#text-measuring); `draw` is one row's cell in `v.columnX` over
// `v.columnWidth`, its batches this module's, handed over per cell (the composer calls it per visible row). An unknown
// type is drawn as `text` (layout and prepare alike).

import { makePaletteFunction } from '@molstar/seq-core/format';
import type { Annotation, Field, SeqSet } from '@molstar/seq-core/model';
import {
  type ColumnRenderer,
  type DrawView,
  fitIn,
  fitted,
  RectBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
} from '../draw/api.js';

/** A column prepared: per alignment row its annotation row (-1: none); per annotation row its text, number and color. */
export interface ColumnState {
  readonly rows: Int32Array;
  readonly text: readonly string[];
  readonly value: Float64Array;
  readonly color: Int32Array;
}
/** A decimal number as text (sign, digits, a point, an exponent), surrounding spaces allowed. */
const DECIMAL = /^\s*[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?\s*$/i;
const fmt = (v: number) => (Number.isInteger(v) ? String(v) : v.toPrecision(3));

/** A column's values formatted once; categories colored by the categorical palette. */
function prepareColumn(kind: string) {
  return (field: Field | null, rows: Int32Array, ann: Annotation | null, set: SeqSet): ColumnState => {
    const n = ann?.count ?? 0;
    const value = new Float64Array(n).fill(Number.NaN);
    const text: string[] = [];
    for (let r = 0; r < n; r++) {
      let s = '';
      if (field?.kind === 'number') value[r] = field.data[r];
      if (field?.kind === 'string' && field.data[r] >= 0) s = set.strings.get(field.data[r]);
      if (field?.kind === 'json' && field.data[r] != null) s = String(field.data[r]);
      // A decimal written out: a blank, `0x1A` or `Infinity` is text, never 0 or a number (GM-07).
      if (field?.kind !== 'number' && DECIMAL.test(s)) value[r] = Number(s);
      if (kind === 'number' || kind === 'fraction')
        s = Number.isNaN(value[r]) ? '' : kind === 'fraction' ? `${Math.round(value[r] * 100)}%` : fmt(value[r]);
      else if (field?.kind === 'number' && !Number.isNaN(value[r])) s = fmt(value[r]);
      text.push(s);
    }
    const color = new Int32Array(n).fill(-1);
    if (kind === 'category') {
      const f = makePaletteFunction(
        { kind: 'categorical', colors: 'Category10' },
        text.filter((t) => t),
      );
      for (let r = 0; r < n; r++) color[r] = text[r] ? (f(text[r]) ?? -1) : -1;
    }
    return { rows, text, value, color };
  };
}

const texts = new TextBatch(1);
const bars = new RectBatch(4);
/** A cell's text box, CSS px (a double passed to a call is boxed): x, width. */
const at = new Float64Array(2);
const fit = new Float64Array(2);

/** `text` cut with an ellipsis to `at[1]` px, placed at `at[0]` aligned `align`, vertically centered in the row. */
function cellText(v: DrawView, text: string, align: number, ink: number): void {
  if (!text || at[1] < 8) return;
  const tm = v.text;
  const px = (v.height * 0.62) | 0;
  const font = tm.font('ui', px < 8 ? 8 : px > 12 ? 12 : px);
  const label = tm.label(text, font);
  fit[0] = at[1];
  fit[1] = tm.label('…', font).width;
  const n = fitIn(label, fit);
  if (n === 0) return;
  texts.reset();
  texts.style(font, align, TextBaseline.middle);
  texts.begin(ink);
  const i = texts.add(fitted(label, n));
  texts.x[i] = Math.round(at[0] * v.dpr);
  texts.y[i] = Math.round((v.top + v.height / 2) * v.dpr) + 0.5;
  v.backend.text(texts);
}
const cell = (s: ColumnState, row: number) => (s && row >= 0 && row < s.rows.length ? s.rows[row] : -1);

export const NumberColumn: ColumnRenderer<ColumnState> = {
  name: 'number',
  width: 64,
  prepare: prepareColumn('number'),
  draw(v, s, row) {
    const r = cell(s, row);
    at[0] = v.columnX + v.columnWidth - 6;
    at[1] = v.columnWidth - 10;
    if (r >= 0) cellText(v, s.text[r], TextAlign.right, SkinToken.ink);
  },
};
export const FractionColumn: ColumnRenderer<ColumnState> = {
  name: 'fraction',
  width: 84,
  prepare: prepareColumn('fraction'),
  draw(v, s, row) {
    const r = cell(s, row);
    if (r < 0 || Number.isNaN(s.value[r])) return;
    const d = v.dpr;
    const x = Math.round((v.columnX + 4) * d) | 0;
    const y = Math.round((v.top + v.height / 2 - 3) * d) | 0;
    const w = Math.round((v.columnWidth - 44) * d) | 0;
    const h = Math.max(1, Math.round(6 * d)) | 0;
    const f = s.value[r];
    bars.reset();
    bars.begin(SkinToken.inkFaint);
    bars.push(x, y, w, h);
    bars.begin(SkinToken.accent);
    bars.push(x, y, Math.round(w * (f < 0 ? 0 : f > 1 ? 1 : f)) | 0, h);
    v.backend.rects(bars);
    at[0] = v.columnX + v.columnWidth - 4;
    at[1] = 34;
    cellText(v, s.text[r], TextAlign.right, SkinToken.inkQuiet);
  },
};
export const TextColumn: ColumnRenderer<ColumnState> = {
  name: 'text',
  width: 110,
  prepare: prepareColumn('text'),
  draw(v, s, row) {
    const r = cell(s, row);
    at[0] = v.columnX + 6;
    at[1] = v.columnWidth - 10;
    if (r >= 0) cellText(v, s.text[r], TextAlign.left, SkinToken.ink);
  },
};
export const CategoryColumn: ColumnRenderer<ColumnState> = {
  name: 'category',
  width: 100,
  prepare: prepareColumn('category'),
  draw(v, s, row) {
    const r = cell(s, row);
    if (r < 0 || !s.text[r]) return;
    const d = v.dpr;
    const sw = v.height < 16 ? v.height - 6 : 10;
    bars.reset();
    bars.begin(s.color[r] >= 0 ? s.color[r] : SkinToken.inkFaint);
    bars.push(
      Math.round((v.columnX + 6) * d) | 0,
      Math.round((v.top + (v.height - sw) / 2) * d) | 0,
      Math.max(1, Math.round(sw * d)) | 0,
      Math.max(1, Math.round(sw * d)) | 0,
    );
    v.backend.rects(bars);
    at[0] = v.columnX + 10 + sw;
    at[1] = v.columnWidth - 16 - sw;
    cellText(v, s.text[r], TextAlign.left, SkinToken.ink);
  },
};
