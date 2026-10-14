/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `emphasis` renderer (seqspec:viewer#extension-points): the covered letters (or an alignment row's cells) of its
// own track, bold in the theme's color or the skin's accent over a cleared cell with an underline, as the track shows
// them (its letters layer's `letter_map`); below the letters LOD its cells in that color, or past 64 items in view
// density per device column as points (each item once at its center, stacked by color: the shared `drawDensity`, with
// its hysteresis). Placed per shown run and clamped to the view or the wrapped line; a ligand's label is kept (an
// underline only) and atoms show their element symbols.
import {
  AxisGeometry,
  CHARS,
  DensityKind,
  type DrawView,
  densityMode,
  drawDensity,
  type EmphasisState,
  elementSymbol,
  type Features,
  itemsInView,
  kernelEntry,
  LETTER_LINE,
  LetterForm,
  Lod,
  letterFontPx,
  letterLines,
  letterRoom,
  NO_COLOR,
  placeRun,
  RectBatch,
  SkinToken,
  segWindow,
  shownHull,
  TextAlign,
  TextBaseline,
  TextBatch,
  ThemeLocation,
  textShown,
  visibleRange,
} from "../draw/api.ts";

const ground = new RectBatch(64);
const marks = new RectBatch();
const text = new TextBatch(64);
const range = new Int32Array(2);
const hull = new Int32Array(2);
const win = new Int32Array(2);
const piece = new Float64Array(16);
let shown = new Int32Array(8);
/** The baselines that center a capital and a lower-case letter on their cell, as the letters draw them. */
const baselines = new Float64Array(2);
const symbol = elementSymbol;

/** Feature k's color: its theme's at its location, else the accent. */
function colorOf(v: DrawView, f: Features, k: number): number {
  const t = v.color;
  if (!t || !v.axis) return SkinToken.accent;
  const loc = ThemeLocation.at(v.loc, v.axis, v.layer, f.block[k], f.start[k]);
  loc.annotation = f.annotation[k];
  loc.row = f.row[k];
  loc.letter = f.code[k] > 0 ? f.code[k] : -1;
  const c = t.value(loc);
  return c === NO_COLOR ? SkinToken.accent : c;
}

function draw(v: DrawView, st: EmphasisState) {
  const ax = v.axis;
  if (!ax) return;
  const f = st.f;
  const dpr = v.dpr;
  const letters = v.geometry === AxisGeometry.wrapped || v.lod === Lod.letters;
  const pad = v.height > 10 ? 2 : 0;
  const y0 = Math.round((v.top + pad) * dpr) | 0;
  const y1 = Math.round((v.top + v.height - pad) * dpr) | 0;
  const under = Math.max(1, Math.round(2 * dpr));
  ground.reset();
  marks.reset();
  text.reset();
  ThemeLocation.reset(v.loc);
  // Past 64 items in view, density; items again below 0.8 of it (one mode per layer per frame).
  if (!letters && densityMode(v, itemsInView(v, f, 65), 64)) return density(v, st, y0, y1);
  for (let b = 0; b < st.letters.codes.length; b++) {
    visibleRange(v, b, range);
    if (range[0] >= range[1]) continue;
    shownHull(v, b, range[0], range[1], hull);
    // On an alignment row (a rows track): that row's letters, where it is held.
    const aln = v.dataRow >= 0 ? ax.set.alignments[st.letters.alignment[b]] : undefined;
    const slot = aln && !aln.lazy ? aln.rows.slotOf(v.dataRow as never) : -1;
    const codes = aln ? (slot >= 0 ? aln.rows.codes : null) : st.letters.codes[b];
    const off = aln ? slot * aln.columnCount : -ax.base[b];
    const form = st.letters.form[b];
    // The letters its track shows (`letter_map`, seqspec:format#layer-types); atoms never mapped.
    const map = form === LetterForm.residues ? st.letters.map : null;
    // Text where the letters (or the row) under it show text (`textShown`), else its cells: never bold letters on a
    // row that shows none (VL-06).
    let rowText = letters;
    if (rowText) {
      letterRoom[0] = v.height;
      const px = v.geometry === AxisGeometry.wrapped ? v.fontPx : letterFontPx(v, b);
      const atoms = form === LetterForm.atoms;
      const plain = v.text.font("mono", px, atoms ? 600 : 400);
      letterRoom[1] = aln ? 0 : atoms ? 2 * v.text.cell(plain, null) + 1 : v.text.cell(plain, st.letters.alphabet[b]);
      rowText = textShown(v, b);
    }
    const shows = rowText && codes !== null;
    if (shows) {
      const font = v.text.font("mono", v.geometry === AxisGeometry.wrapped ? v.fontPx : letterFontPx(v, b), 600);
      if (text.n > 0 && text.font !== font) {
        v.backend.text(text);
        text.reset();
      }
      text.style(font, TextAlign.center, TextBaseline.alphabetic);
      letterLines(v, font, y0, y1, baselines);
    }
    for (let lane = 0; lane < f.lanes; lane++) {
      const g = b * f.lanes + lane;
      segWindow(f.start, f.reach, f.seg[g], f.seg[g + 1], hull[0], hull[1], win);
      for (let k = win[0]; k < win[1]; k++) {
        if (Math.max(f.stop[k], f.start[k] + 1) <= hull[0] || (f.on[k] >= 0 && f.on[k] !== v.row)) continue;
        const color = colorOf(v, f, k);
        marks.begin(color);
        if (!shows || !codes) {
          // A row showing text this layer cannot read (a streamed block's or alignment's letters): the underline
          // alone, never a cell over the letters (VL-10); a row without text: its cells.
          const top = rowText ? y1 - under : y0;
          const n = placeRun(v, b, f.start[k], Math.max(f.stop[k], f.start[k] + 1), piece);
          for (let q = 0; q < n; q++)
            marks.push(piece[2 * q], top, Math.max(1, piece[2 * q + 1] - piece[2 * q]), y1 - top);
          continue;
        }
        let m = v.runToShown(b, f.start[k], Math.max(f.stop[k], f.start[k] + 1), shown);
        if (2 * m > shown.length) {
          shown = new Int32Array(4 * m);
          m = v.runToShown(b, f.start[k], Math.max(f.stop[k], f.start[k] + 1), shown);
        }
        text.begin(color);
        for (let r = 0; r < m; r++)
          for (let q = Math.max(range[0], shown[2 * r]); q < Math.min(range[1], shown[2 * r + 1]); q++) {
            const x = v.xd(b, q);
            const w = v.xe(b, q) - x;
            marks.push(x, y1 - under, w, under);
            if (form === LetterForm.ligand) continue;
            ground.begin(SkinToken.surface);
            ground.push(x, y0, w, y1 - y0);
            const c0 = codes[off + ax.indexOf(b, q as never)];
            const c = map !== null ? map[c0 & 127] : c0;
            const t = text.add(form === LetterForm.atoms ? symbol(c) : CHARS[c & 127]);
            text.x[t] = x + w / 2;
            text.y[t] = baselines[form === LetterForm.atoms ? 0 : LETTER_LINE[c & 127]];
          }
      }
    }
  }
  if (ground.n > 0) v.backend.rects(ground);
  if (marks.n > 0) v.backend.rects(marks);
  if (text.n > 0) v.backend.text(text);
}

/** Density as points: the shared band, each item once at its center, stacked by color. */
function density(v: DrawView, st: EmphasisState, y0: number, y1: number) {
  marks.reset();
  drawDensity(v, st.f, st, true, DensityKind.stacked, SkinToken.accent, marks, y0, y1);
  if (marks.n > 0) v.backend.rects(marks);
}

export const EmphasisRenderer = kernelEntry<EmphasisState>("emphasis", 22, {
  geometries: ["linear", "wrapped"],
  stable: true,
  draw,
  tiny: (v, s) => density(v, s, Math.round(v.top * v.dpr) | 0, Math.round((v.top + v.height) * v.dpr) | 0),
});
