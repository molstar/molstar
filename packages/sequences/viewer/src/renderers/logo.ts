/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `logo` renderer (seqspec:viewer#extension-points; seqspec:sources#levels-and-the-zoom-gate): each position's
// letters stacked by bits or probability, lowest first, each letter's height its share of the band (a column's stack is
// walked once with a running `below` offset, never once per letter). A letter is a character stretched over
// its height and its cell (the text batch's per-item scale), as readers expect of a sequence logo; below two device px
// a position, and on the tiny path, one bar per device column of the tallest stack there, in its top letter's color,
// from a per-block pyramid of the tallest stack per bin, so a frame costs its device columns, not its positions.
import {
  type Axis,
  AxisGeometry,
  CHARS,
  columnBins,
  columnEnd,
  columnLevel,
  type DrawView,
  kernelEntry,
  type MatrixState,
  NO_COLOR,
  PYRAMID_BASE,
  RectBatch,
  SkinToken,
  TextAlign,
  TextBaseline,
  TextBatch,
  ThemeLocation,
  visibleRange,
} from '../draw/api.js';

const bars = new RectBatch(256);
const letters = new TextBatch(256);
/** The letters' font size (CSS px) and its capitals' height (about 0.72 of it), which a letter's scale stretches. */
const PX = 20;
const CAP = 0.72 * PX;
let colors = new Int32Array(32);
const range = new Int32Array(2);
const bins = new Int32Array(2);
/** The band's bottom and height (device px): doubles, so they reach the walks through a typed slot, never a call. */
const frame = new Float64Array(2);

/** A block's stacks for the thin path: per position its height (0 without a row), per level the tallest per bin. */
interface Tallest {
  total: Float32Array;
  /** Level `PYRAMID_BASE + j`: per bin of 2^level positions, the position of its tallest stack (ties: the first). */
  levels: Int32Array[];
}
const tallest = new WeakMap<MatrixState, (Tallest | null)[]>();
const tallestOf = (s: MatrixState) => tallest.get(s) ?? newTallest(s);
/** Edge-rate, apart from the draw: built on the first thin draw of a prepared state. */
function newTallest(s: MatrixState): (Tallest | null)[] {
  const S = s.states;
  const x = s.heights.map((h, b) => {
    const rows = s.rows[b];
    if (!h || !rows) return null;
    const n = rows.length;
    const total = new Float32Array(n);
    for (let i = 0; i < n; i++) for (let q = 0; rows[i] >= 0 && q < S; q++) total[i] += h[i * S + q];
    const levels: Int32Array[] = [];
    for (let size = 1 << PYRAMID_BASE, below: Int32Array | null = null; size < 2 * n; size *= 2) {
      const at = new Int32Array(Math.ceil(n / size));
      for (let j = 0; j < at.length; j++) {
        // From the level below (two halves), or at the base from its positions.
        let best = below ? below[2 * j] : j * size;
        const last = below ? Math.min(2 * j + 2, below.length) : Math.min((j + 1) * size, n);
        for (let c = below ? 2 * j + 1 : j * size + 1; c < last; c++) {
          const i = below ? below[c] : c;
          if (total[i] > total[best]) best = i;
        }
        at[j] = best;
      }
      levels.push(at);
      below = at;
    }
    return { total, levels };
  });
  tallest.set(s, x);
  return x;
}

/** Block b's letters at shown [range[0], range[1]), each stack walked once from the bottom, stretched to its cell. */
function stacks(v: DrawView, s: MatrixState, b: number) {
  const ax = v.axis as Axis;
  const S = s.states;
  const h = s.heights[b] as Float32Array;
  const o = s.order[b] as Uint8Array;
  const rows = s.rows[b] as Int32Array;
  const adv = v.text.advances(letters.font);
  const dpr = v.dpr;
  for (let k = range[0]; k < range[1]; k++) {
    const i = ax.indexOf(b, k as never) - ax.base[b];
    if (rows[i] < 0) continue;
    const x0 = v.xd(b, k);
    const w = v.xe(b, k) - x0;
    // Half a CSS px clear on each side, as the old stack's letters (a bar's width where the cell is narrower).
    const lw = Math.max(1, w - dpr);
    let below = 0;
    for (let r = 0; r < S; r++) {
      const q = o[i * S + r];
      const hq = h[i * S + q] * frame[1];
      if (hq >= dpr) {
        const code = s.codes[q] & 127;
        letters.begin(colors[q]);
        const j = letters.add(CHARS[code]);
        letters.x[j] = x0 + (w - lw) / 2;
        letters.y[j] = frame[0] - below;
        letters.sx[j] = lw / ((adv[code] || 12) * dpr);
        letters.sy[j] = hq / (CAP * dpr);
      }
      below += hq;
    }
  }
}

/** Block b's shown [range[0], range[1]) as one bar per device column: the tallest stack there, its top letter's color. */
function columnsOf(v: DrawView, s: MatrixState, b: number) {
  const ax = v.axis as Axis;
  const S = s.states;
  const t = tallestOf(s)[b] as Tallest;
  const total = t.total;
  const o = s.order[b] as Uint8Array;
  const level = columnLevel(v, b, t.levels.length);
  const bottom = frame[0] | 0;
  for (let k = range[0]; k < range[1]; ) {
    const k1 = columnEnd(v, b, k);
    let best = -1;
    if (level >= 0 && columnBins(v, b, k, k1, level, bins)) {
      const at = t.levels[level - PYRAMID_BASE];
      for (let j = bins[0]; j < bins[1]; j++) if (best < 0 || total[at[j]] > total[best]) best = at[j];
    } else
      for (let j = k; j < k1; j++) {
        const i = ax.indexOf(b, j as never) - ax.base[b];
        if (best < 0 || total[i] > total[best]) best = i;
      }
    if (best >= 0 && total[best] > 0) {
      bars.begin(colors[o[best * S + S - 1]]);
      const x = v.xd(b, k);
      const hb = Math.max(1, Math.round(total[best] * frame[1]));
      bars.push(x, bottom - hb, Math.max(1, v.xd(b, k1) - x), hb);
    }
    k = k1;
  }
}

function draw(v: DrawView, s: MatrixState, tiny: boolean) {
  const ax = v.axis;
  bars.reset();
  letters.reset();
  if (!ax || v.geometry !== AxisGeometry.linear || s.states === 0) return;
  if (colors.length < s.states) colors = new Int32Array(2 * s.states);
  letters.style(v.text.font('mono', PX, 600), TextAlign.left, TextBaseline.alphabetic);
  const pad = v.height >= 8 ? 2 : 0;
  frame[0] = Math.round((v.top + v.height - pad) * v.dpr);
  frame[1] = (v.height - 2 * pad) * v.dpr;
  for (let b = 0; frame[1] >= 1 && b < ax.count; b++) {
    visibleRange(v, b, range);
    if (!s.heights[b] || !s.rows[b] || range[0] >= range[1]) continue;
    // Each state's color at this block (a letter theme: by its letter).
    ThemeLocation.reset(v.loc);
    ThemeLocation.at(v.loc, ax, v.layer, b, ax.base[b]);
    for (let q = 0; q < s.states; q++) {
      v.loc.letter = s.codes[q];
      const c = v.color ? v.color.value(v.loc) : NO_COLOR;
      colors[q] = c >= 0 ? c : SkinToken.ink;
    }
    if (tiny || v.cell * v.span[b] * v.dpr < 2) columnsOf(v, s, b);
    else stacks(v, s, b);
  }
  if (bars.n > 0) v.backend.rects(bars);
  if (letters.n > 0) v.backend.text(letters);
}

export const LogoRenderer = kernelEntry<MatrixState>('logo', 48, {
  minHeight: 8,
  draw: (v, s) => draw(v, s, false),
  tiny: (v, s) => draw(v, s, true),
});
