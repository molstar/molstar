/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `letters` renderer (seqspec:viewer#extension-points): residues, then cells, thin, a bar; a ligand as one labeled
// cell; `row`: an alignment row's letters; `of`: other polymers' residues through each placement's path, one letter or
// three by the room (with a margin each side of the part) and `residue_form`, each over its larger part. Cells by the color theme, reduced by its kind below
// a device px (seqspec:sources#levels-and-the-zoom-gate); letters where the block's widest letter fits, in the ink that
// reads on their cell, atoms by element symbol; wrapped, `inlineStart` group labels before each group. Every location
// is filled with `ThemeLocation.at`. Batches are this module's, reused by every draw. ST3a: a block over a lazy
// molecule draws its letters from the stream's regions (`StreamView.letters`), copied into module scratch for the shown
// hull, where a device px holds at most `LETTERS_PER_PX` positions; cells thinner than a device px read raw (no index
// over a chromosome), and positions not yet held draw nothing: a pending region's placeholder band, and past the
// letters gate a hatched band with its hint (ST3a-2, `draw/streamed.ts`).
// A row with no cell colors (text only), where its text does not fit, marks its residues: dots, then bars (`DOT_OF_CELL`).

import { CoordinateSystemKind, threeLetter } from '@molstar/seq-core/model';
import {
  type Axis,
  AxisGeometry,
  type AxisView,
  bandOver,
  type CellSource,
  CHARS,
  cellBand,
  cellInk,
  type DrawView,
  drawCells,
  elementSymbol,
  type FontId,
  fitIn,
  fitted,
  GlyphBatch,
  GlyphShape,
  Granularity,
  type Hit,
  HoverKind,
  inlineGroupLabels,
  kernelEntry,
  LETTER_LINE,
  LetterForm,
  type LettersState,
  letterCap,
  letterFontPx,
  letterLines,
  letterRoom,
  lettersGate,
  NO_COLOR,
  type OfLetters,
  OPAQUE,
  ofPartAt,
  PathBatch,
  placeRun,
  RectBatch,
  RegionState,
  SkinToken,
  segWindow,
  shownHull,
  streamedBlock,
  TextAlign,
  TextBaseline,
  TextBatch,
  type TextManagerApi,
  ThemeLocation,
  textShown,
  visibleRange,
  type Skin,
} from '../draw/api.js';

type Params = Readonly<Record<string, unknown>>;

/** Unstable with `of` or a ligand form, whose labels are wider than a cell. */
function stable(s: LettersState): boolean {
  if (s.of) return false;
  for (let b = 0; b < s.form.length; b++) if (s.form[b] === LetterForm.ligand) return false;
  return true;
}

/** `of` another polymer: its residue's location; `row`: the row's residue. */
function hit(v: AxisView, s: LettersState, _p: unknown, out: Hit): boolean {
  const o = s.of;
  if (o && out.block >= 0) {
    const k = out.block + 1 < o.seg.length ? ofPartAt(o, out.block, out.index) : -1;
    if (k < 0 || o.unit[k] < 0) return false;
    out.kind = HoverKind.seq;
    out.seq.unit = o.unit[k] as never;
    out.seq.index = o.source[k] as never;
    out.system.kind = CoordinateSystemKind.unit;
    out.system.id = o.unit[k];
    out.system.row = -1;
    return true;
  }
  if (s.row < 0 || out.block < 0 || s.alignment[out.block] < 0) return false;
  const set = (v.axis as Axis).set;
  const aln = set.alignments[s.alignment[out.block]];
  const index = aln ? aln.rows.colToPos(s.row as never, out.index as never) : -1;
  if (!aln || index < 0) return false;
  out.kind = HoverKind.row;
  out.row.set = set;
  out.row.alignment = aln.id;
  out.row.row = s.row as never;
  out.row.index = index as never;
  out.system.kind = CoordinateSystemKind.row;
  out.system.id = aln.id;
  out.system.row = s.row;
  return true;
}

const rects = new RectBatch();
/** Hairlines where a codon starts, over the cells (a batch of their own: a backend groups a batch by color). */
const lines = new RectBatch(64);
const letters = new TextBatch();
const sups = new TextBatch(32);
const pills = new TextBatch(16);
/** Past the letters gate: the hatch and its hint. */
const hatch = new PathBatch(256);
const hints = new TextBatch(8);
const range = new Int32Array(2);
const hull = new Int32Array(2);
const win = new Int32Array(2);
const piece = new Float64Array(8);
/** The margin (CSS px) a residue's three-letter code keeps from each edge of the part it is centered over. */
const OF_MARGIN = 2;
const fit = new Float64Array(2);
/**
 * Where text does not fit and a position has no cell color (a text-only row: `of` with no cell color, a letters row
 * with none), its presence is marked: a dot in its cell (a diameter of `DOT_OF_CELL` of the cell, at most `DOT_MAX` CSS
 * px and the row's share), else, under `DOT_MIN` device px, a thin bar over each run of present positions, merged per
 * device column. One form for the whole row.
 */
const DOT_OF_CELL = 0.4;
const DOT_MAX = 9;
const DOT_MIN = 3;
/** The marks' strength (of 255): the ink at half, quieter than the letters they stand for, in every skin. */
const MARK_ALPHA = 128;
const dots = new GlyphBatch(128);
const bars = new RectBatch(64);
const MarkForm = { none: 0, dots: 1, bars: 2 } as const;
const MARK_NONE = MarkForm.none;
const MARK_DOTS = MarkForm.dots;
const MARK_BARS = MarkForm.bars;
const CIRCLE = GlyphShape.circle;
const INK = SkinToken.ink;
/** The mark geometry (device px): [0] the dot's diameter, from the cell width in [0] on entry; [1] the row's center; [2] the bar's top, [3] its height. */
const mk = new Float64Array(4);
/** The open bar's edges (device px); its end before its start: none open. */
const open = new Float64Array(2);

/** `mk`'s geometry for marks over a cell `mk[0]` CSS px wide in the band [y0, y1) (device px); the form it gives. */
function markForm(v: DrawView, y0: number, y1: number): number {
  const dpr = v.dpr;
  // Integer and branch arithmetic: no Math builtin (each boxed a double under the test runner's module wrapper).
  let d = DOT_OF_CELL * mk[0];
  if (d > DOT_MAX) d = DOT_MAX;
  const cap = (v.height * 0.6) | 0;
  if (d > cap) d = cap;
  mk[0] = d * dpr;
  mk[1] = (y0 + y1) / 2;
  let h = ((y1 - y0) / 4 + 0.5) | 0;
  if (h > 3 * dpr) h = (3 * dpr + 0.5) | 0;
  if (h < 1) h = 1;
  mk[3] = h;
  mk[2] = ((y0 + y1 - h) / 2 + 0.5) | 0;
  return mk[0] >= DOT_MIN ? MARK_DOTS : MARK_BARS;
}
/** A dot centered over the device columns [x0, x1). */
function pushDot(x0: number, x1: number): void {
  const i = dots.add(CIRCLE);
  dots.x[i] = (x0 + x1) / 2;
  dots.y[i] = mk[1];
  dots.radius[i] = mk[0] / 2;
}
/** Extends the open bar over [x0, x1), or closes it and opens one there (a device px of rounding bridged). */
function pushBar(x0: number, x1: number): void {
  if (open[1] >= open[0] && x0 <= open[1] + 1) {
    if (x1 > open[1]) open[1] = x1;
    return;
  }
  closeBar();
  open[0] = x0;
  open[1] = x1;
}
function closeBar(): void {
  if (open[1] >= open[0]) bars.push(open[0], mk[2], Math.max(1, open[1] - open[0]), mk[3]);
  open[0] = 0;
  open[1] = -1;
}
/** The baselines that center a capital and a lower-case letter on their cell (`letterLines`). */
const baselines = new Float64Array(2);
/** Three-letter codes by one-letter code, made on first use; element symbols are `elementSymbol`'s. */
const THREE: string[] = [];
const symbol = elementSymbol;
const three = (c: number) => {
  THREE[c] ??= threeLetter(c);
  return THREE[c];
};

/** A streamed block's codes over its shown hull (`StreamView.letters`), grown at edge rate. */
let held = new Uint8Array(4096);
/**
 * A block's own letters as cells: the color theme at block index i, keyed by letter where it colors by letter alone.
 * `codes` and `base` are the block's being drawn (a streamed one's: `held` from its hull's start; `raw`: no index).
 */
/** The skin's `unobserved` color as a number, parsed when the skin changes (a draw reads it per position). */
let grayOf: Skin | null = null;
let gray = NO_COLOR;
function unobservedColor(skin: Skin): number {
  if (skin !== grayOf) {
    grayOf = skin;
    const m = /^rgba?\(\s*(\d+)[\s,]+(\d+)[\s,]+(\d+)/.exec(skin.unobserved);
    const hex = /^#([0-9a-f]{6})$/i.exec(skin.unobserved);
    gray = m ? (Number(m[1]) << 16) | (Number(m[2]) << 8) | Number(m[3]) : hex ? Number.parseInt(hex[1], 16) : NO_COLOR;
  }
  return gray;
}
/** Whether block b's position i is a unit position the structure does not model (the unit's `observed` ranges). */
function unobservedAt(ax: Axis, b: number, i: number): boolean {
  const u = ax.unit[b];
  const obs = u >= 0 ? ax.set.units.observed[u] : null;
  if (!obs) return false;
  for (let k = 0; k < obs.length; k += 2) if (i >= obs[k] && i < obs[k + 1]) return false;
  return true;
}

const own: CellSource & {
  v: DrawView | null;
  s: LettersState | null;
  codes: Uint8Array | null;
  base: number;
  raw: boolean;
} = {
  v: null,
  s: null,
  codes: null,
  base: 0,
  raw: false,
  color(b, i) {
    const v = own.v as DrawView;
    const s = own.s as LettersState;
    const codes = own.codes;
    if (!v.color || !codes || !v.axis) return NO_COLOR;
    if (v.props.unobservedGray && unobservedAt(v.axis, b, i)) return unobservedColor(v.skin);
    // The theme reads the letter shown (`letter_map`), not the stored one; atoms are numbers, never mapped.
    const c = codes[i - own.base];
    if (c === 0) return NO_COLOR;
    const loc = ThemeLocation.at(v.loc, v.axis, v.layer, b, i);
    loc.row = s.row;
    loc.letter = s.map !== null && s.form[b] === LetterForm.residues ? s.map[c & 127] : c;
    return v.color.value(loc);
  },
  alpha() {
    const t = (own.v as DrawView).color;
    return t?.alpha ? t.alpha((own.v as DrawView).loc) : OPAQUE;
  },
  key(_b, i) {
    const v = own.v as DrawView;
    const g = v.color ? v.color.granularity : -1;
    const codes = own.codes;
    if (g === Granularity.letter && codes && v.axis) return codes[i - own.base];
    return g === Granularity.uniform || g === Granularity.unit ? 0 : -1;
  },
  meanValues: () => null,
  meanPyramid: () => null,
};
/** Another polymer's residues as cells: each residue's color over its part's cells, keyed by part. */
const ofCells: CellSource & { v: DrawView | null; o: OfLetters | null; map: Uint8Array | null } = {
  v: null,
  o: null,
  map: null,
  color(b, i) {
    const e = ofPartAt(ofCells.o as OfLetters, b, i);
    return e < 0 ? NO_COLOR : ofColor(ofCells.v as DrawView, ofCells.o as OfLetters, b, e, ofCells.map);
  },
  alpha() {
    const t = (ofCells.v as DrawView).color;
    return t?.alpha ? t.alpha((ofCells.v as DrawView).loc) : OPAQUE;
  },
  key: (b, i) => ofPartAt(ofCells.o as OfLetters, b, i) + 1,
  // Its parts' hull on the block: a protein's codons on its chromosome, never the chromosome (GM-02).
  extent(b, out) {
    const o = ofCells.o as OfLetters;
    const lo = o.seg[b];
    const hi = o.seg[b + 1];
    out[0] = hi > lo ? o.start[lo] : 0;
    out[1] = hi > lo ? o.reach[hi - 1] : 0;
  },
  meanValues: () => null,
  meanPyramid: () => null,
};
/** Part e's residue in the theme (its polymer's unit, its own index and letter). */
function ofColor(v: DrawView, o: OfLetters, b: number, e: number, map: Uint8Array | null): number {
  if (!v.color || !v.axis) return NO_COLOR;
  const loc = ThemeLocation.at(v.loc, v.axis, v.layer, b, o.start[e]);
  loc.unit = o.unit[e];
  loc.index = o.source[e];
  loc.letter = map !== null ? map[o.code[e] & 127] : o.code[e];
  return v.color.value(loc);
}

/**
 * `textPx`'s cell width, row height and cap (`letterCap`), CSS px: module scratch, so no double is boxed across the
 * call.
 */
const room = new Float64Array(3);
/**
 * The font px of `chars` characters to a cell `room[0]` CSS px wide in a row `room[1]` tall, the letters' size (never
 * larger than the row allows, nor than the cap `room[2]`); every residue text of a row is drawn at it.
 */
function textPx(chars: number): number {
  return Math.max(6, Math.min(room[2], Math.floor(Math.min((room[0] / chars) * 1.25, room[1] * 0.72)))) | 0;
}
/** A font's generic cell width (CSS px), `monoCell`'s result: module scratch, as `room`. */
const cellW = new Float64Array(1);
/**
 * `tm.cell(font, null)` into `cellW[0]`. `drawOf` is too large for TurboFan to inline `cell` at its own call sites, so
 * the width came back boxed (16 bytes a call, two a frame on case 14's codon row); this helper is small enough to be
 * inlined with `cell` inside it, and no double crosses a call.
 */
function monoCell(tm: TextManagerApi, font: FontId): void {
  cellW[0] = tm.cell(font, null);
}
/** Block b's letters' font px. */
function letterPx(v: DrawView, b: number): number {
  if (v.geometry === AxisGeometry.wrapped) return v.fontPx;
  letterRoom[0] = v.height;
  return letterFontPx(v, b);
}

/** Hands the letters batch over when the next block's font differs (one font per batch). */
function letterFont(v: DrawView, font: number) {
  if (letters.n > 0 && letters.font !== font) {
    v.backend.text(letters);
    letters.reset();
  }
  letters.style(font, TextAlign.center, TextBaseline.alphabetic);
}

function draw(v: DrawView, s: LettersState, p: unknown) {
  const ax = v.axis;
  if (!ax) return;
  rects.reset();
  lines.reset();
  dots.reset();
  bars.reset();
  bars.begin(INK, MARK_ALPHA);
  dots.begin(INK, MARK_ALPHA);
  open[1] = -1;
  letters.reset();
  sups.reset();
  pills.reset();
  hatch.reset();
  hints.reset();
  ThemeLocation.reset(v.loc);
  try {
    if (s.of) drawOf(v, s, s.of, p as Params);
    else for (let b = 0; b < s.codes.length; b++) block(v, s, b);
    if (rects.n > 0) v.backend.rects(rects);
    if (lines.n > 0) v.backend.rects(lines);
    if (bars.n > 0) v.backend.rects(bars);
    if (dots.n > 0) {
      dots.font = v.text.font('ui', 12, 600);
      v.backend.glyphs(dots);
    }
    if (letters.n > 0) v.backend.text(letters);
    if (sups.n > 0) v.backend.text(sups);
    if (pills.n > 0) v.backend.text(pills);
    if (hatch.n > 0) v.backend.paths(hatch);
    if (hints.n > 0) v.backend.text(hints);
  } finally {
    own.v = own.s = own.codes = ofCells.v = ofCells.o = ofCells.map = null;
  }
}

/**
 * A streamed block's codes over the hull of shown [lo, hi) into `held`: the hull's first index, or -1 where nothing is
 * held. Past the letters gate the hatched band (nothing asked there); a run not held whose region is pending, its
 * placeholder band, between device rows y0 and y1.
 */
function streamed(v: DrawView, b: number, lo: number, hi: number, y0: number, y1: number): number {
  const st = v.stream;
  if (!st?.letters || !streamedBlock(v, b) || lettersGate(v, b, lo, hi, y0, y1, hatch, hints)) return -1;
  shownHull(v, b, lo, hi, hull);
  const n = hull[1] - hull[0];
  if (held.length < n) held = new Uint8Array(Math.max(n, 2 * held.length));
  const got = st.letters(b, hull[0], hull[1], held);
  for (let a = 0; got < n && a < n; ) {
    let e = a;
    while (e < n && held[e] === 0) e++;
    if (e > a && st.stateAt?.(b, hull[0] + a) === RegionState.pending)
      bandOver(v, b, hull[0] + a, hull[0] + e, y0, y1, SkinToken.pending, rects);
    a = e + 1;
  }
  return got > 0 ? hull[0] : -1;
}

/** Positions past which a text-only block's presence is one bar over its shown range (gaps are not looked for). */
const MARK_SCAN = 1 << 15;
/**
 * A block's positions marked where text does not fit and the row has no cell color (see `DOT_OF_CELL`): a dot in each
 * present position's cell, or a bar over each run of them; a gap column or an unheld position is absent, a break.
 * `own` holds the block.
 */
function blockMarks(v: DrawView, b: number, lo: number, hi: number, base: number, atoms: boolean) {
  if (v.geometry === AxisGeometry.wrapped || v.height < 4) return;
  const ax = v.axis as Axis;
  const codes = own.codes as Uint8Array;
  let k = lo;
  while (k < hi) {
    const c = codes[ax.indexOf(b, k as never) - base];
    if (c !== 0 && (atoms || (c !== 45 && c !== 46))) break;
    k++;
  }
  if (k === hi || (v.color !== null && own.color(b, ax.indexOf(b, k as never)) >= 0)) return;
  const pad = v.height > 10 ? 2 : 0;
  mk[0] = v.cell * v.span[b];
  const form = markForm(v, Math.round((v.top + pad) * v.dpr) | 0, Math.round((v.top + v.height - pad) * v.dpr) | 0);
  // Decided over the whole view (`viewLo`), so a pan's strip draws what a full frame draws there (VL-07).
  if (v.viewHi[b] - v.viewLo[b] > MARK_SCAN) {
    pushBar(v.xd(b, lo), v.xe(b, hi - 1));
    return closeBar();
  }
  for (; k < hi; k++) {
    const i = ax.indexOf(b, k as never);
    const c = codes[i - base];
    if (c === 0 || (!atoms && (c === 45 || c === 46)) || own.color(b, i) >= 0) continue;
    if (form === MARK_DOTS) pushDot(v.xd(b, k), v.xe(b, k));
    else pushBar(v.xd(b, k), v.xe(b, k));
  }
  if (form === MARK_BARS) closeBar();
}

function block(v: DrawView, s: LettersState, b: number) {
  const ax = v.axis as Axis;
  visibleRange(v, b, range);
  const lo = range[0];
  const hi = range[1];
  let codes = s.codes[b];
  let base = ax.base[b];
  own.raw = false;
  if (!codes && lo < hi) {
    const pad = v.height > 10 ? 2 : 0;
    base = streamed(
      v,
      b,
      lo,
      hi,
      Math.round((v.top + pad) * v.dpr) | 0,
      Math.round((v.top + v.height - pad) * v.dpr) | 0,
    );
    codes = base < 0 ? null : held;
    own.raw = true;
  }
  if (!codes || lo >= hi) return;
  const form = s.form[b];
  if (form === LetterForm.ligand) return pill(v, s, b, lo, hi);
  const wrapped = v.geometry === AxisGeometry.wrapped;
  const w = v.cell * v.span[b];
  // Residues and atoms alike: their cells in the theme's color at every zoom (a ligand's one position is its pill).
  const cells = !wrapped && v.color !== null;
  // Wrapped text with `wrappedColors`: the letters themselves in the layer's color (there are no cells to put it on).
  const tinted = wrapped && v.color !== null && v.props.wrappedColors;
  own.v = v;
  own.s = s;
  own.codes = codes;
  own.base = base;
  const pad = v.height > 10 ? 2 : 0;
  if (cells) {
    cellBand[0] = Math.round((v.top + pad) * v.dpr);
    cellBand[1] = Math.round((v.top + v.height - pad) * v.dpr) - cellBand[0];
    // `cellGap`'s one device px right of a wide cell (inline: a double passed to a call is boxed below TurboFan).
    cellBand[2] = w > 6 ? 1 : 0;
    drawCells(v, b, lo, hi, own, s, rects);
  }
  const atoms = form === LetterForm.atoms;
  const font = v.text.font('mono', letterPx(v, b), atoms ? 600 : 400);
  // The letters LOD per block: the block's widest letter, or two for an element symbol.
  const need = atoms ? 2 * v.text.cell(font, null) + 1 : v.text.cell(font, s.alphabet[b]);
  letterRoom[0] = v.height;
  letterRoom[1] = need;
  if (!textShown(v, b)) return blockMarks(v, b, lo, hi, base, atoms);
  if (wrapped) {
    sups.style(v.text.font('mono', Math.max(6, Math.round(v.fontPx * 0.62))), TextAlign.right, TextBaseline.middle);
    inlineGroupLabels(v, b, lo, hi, sups);
  }
  letterFont(v, font);
  const ya = Math.round((v.top + pad) * v.dpr) | 0;
  letterLines(v, font, ya, Math.round((v.top + v.height - pad) * v.dpr) | 0, baselines);
  const map = atoms ? null : s.map;
  for (let k = lo; k < hi; k++) {
    const i = ax.indexOf(b, k as never);
    const stored = codes[i - base];
    // A streamed position not yet held.
    if (stored === 0 && !atoms) continue;
    const c = map !== null ? map[stored & 127] : stored;
    // An atom in its theme's color; a residue in the ink that reads on its cell; a gap quiet.
    const fill = atoms || cells || tinted ? own.color(b, i) : NO_COLOR;
    letters.begin(
      atoms
        ? fill >= 0
          ? fill
          : SkinToken.ink
        : c === 45 || c === 46
          ? SkinToken.inkQuiet
          : fill >= 0
            ? tinted
              ? fill
              : cellInk(fill)
            : SkinToken.ink,
    );
    const t = letters.add(atoms ? symbol(c) : CHARS[c & 127]);
    letters.x[t] = (v.xd(b, k) + v.xe(b, k)) / 2;
    letters.y[t] = baselines[atoms ? 0 : LETTER_LINE[c & 127]];
  }
}

/** A ligand in the `position` form: one cell across its span, labeled with its chemical component. */
function pill(v: DrawView, s: LettersState, b: number, lo: number, hi: number) {
  const dpr = v.dpr;
  const d = Math.max(1, Math.round(dpr));
  const y0 = Math.round((v.top + 3) * dpr);
  const y1 = Math.round((v.top + v.height - 3) * dpr);
  const tm = v.text;
  const font = tm.font('ui', Math.max(7, Math.min(12, Math.floor(v.height * 0.5))), 600);
  const name = tm.label(s.names[b] ?? '', font);
  pills.style(font, TextAlign.center, TextBaseline.alphabetic);
  letterLines(v, font, y0, y1, baselines);
  pills.begin(SkinToken.ink);
  for (let k = lo; k < hi; k++) {
    const x0 = v.xd(b, k) + d;
    const x1 = v.xe(b, k) - d;
    if (x1 - x0 < 2 * d) continue;
    rects.begin(SkinToken.inkFaint);
    rects.push(x0, y0, x1 - x0, y1 - y0);
    rects.begin(SkinToken.surfaceAlt);
    rects.push(x0 + d, y0 + d, x1 - x0 - 2 * d, y1 - y0 - 2 * d);
    fit[0] = (x1 - x0) / dpr - 6;
    fit[1] = tm.label('…', font).width;
    const n = fitIn(name, fit);
    if (n === 0) continue;
    const t = pills.add(fitted(name, n));
    pills.x[t] = (x0 + x1) / 2;
    pills.y[t] = baselines[0];
  }
}

/**
 * Another polymer's residues (`of`): each residue's cells in its color (thinner than a device px, reduced as a block's
 * cells), a hairline where a residue over several positions (a codon) starts, its text centered over its larger part's
 * widest piece. The text follows the letters' rule, one size and one form for the whole row: three letters (`Met`)
 * where a residue spans a codon, `residue_form` allows and the narrowest larger part in view holds them with a margin
 * (`OF_MARGIN`) each side, else one letter where it holds the alphabet's widest, else none; never a residue's text
 * left out alone, but for one the view's edge cuts narrower than its text. Placed through `placeRun`, so a
 * wrapped line clamps them.
 */
function drawOf(v: DrawView, s: LettersState, o: OfLetters, p: Params) {
  const tm = v.text;
  const dpr = v.dpr;
  const pad = v.height > 10 ? 2 : 0;
  const y0 = Math.round((v.top + pad) * dpr);
  const y1 = Math.round((v.top + v.height - pad) * dpr);
  const linear = v.geometry === AxisGeometry.linear;
  const form = p.residue_form;
  ofCells.v = v;
  ofCells.o = o;
  ofCells.map = s.map;
  for (let b = 0; b + 1 < o.seg.length; b++) {
    visibleRange(v, b, range);
    if (o.seg[b + 1] === o.seg[b] || range[0] >= range[1]) continue;
    const thin = linear && v.cell * v.span[b] * dpr < 1;
    if (thin) {
      cellBand[0] = y0;
      cellBand[1] = y1 - y0;
      cellBand[2] = 0;
      drawCells(v, b, range[0], range[1], ofCells, s, rects);
    }
    shownHull(v, b, range[0], range[1], hull);
    segWindow(o.start, o.reach, o.seg[b], o.seg[b + 1], hull[0], hull[1], win);
    // The positions a residue spans here (3 over its codon), and the room its text has: the row's one rule.
    // The text is centered over a residue's larger part, so the room is the narrowest larger part's (a codon an
    // intron splits has two cells there, not three), not the widest residue's.
    let span = 1;
    let least = 1 << 30;
    for (let e = win[0]; e < win[1]; e++) {
      if (o.main[e] === 0) continue;
      const n = o.stop[e] - o.start[e];
      if (n > span) span = n;
      if (n < least) least = n;
    }
    if (least === 1 << 30) least = span;
    if (!linear) monoCell(tm, tm.font('mono', v.fontPx));
    room[0] = least * (linear ? v.cell * v.span[b] : cellW[0]);
    room[1] = v.height;
    room[2] = letterCap(v.props);
    let chars = 0;
    let font = 0;
    if (!linear || (v.height >= 10 && room[0] >= v.props.minLetterCell)) {
      if (form === 'three' || (form !== 'one' && span >= 3)) {
        font = tm.font('mono', linear ? textPx(3) : v.fontPx);
        // With a margin each side of the larger part (and the hairline that starts a codon): the whole row or none.
        monoCell(tm, font);
        if (room[0] >= 3 * cellW[0] + 2 * OF_MARGIN + (linear && least > 1 ? 1 : 0)) chars = 3;
      }
      if (chars === 0 && form !== 'three') {
        font = tm.font('mono', linear ? textPx(1) : v.fontPx);
        monoCell(tm, font);
        if (room[0] >= cellW[0]) chars = 1;
      }
    }
    // Text-only row (its first residue has no cell color) where no text fits: a dot or a bar of presence.
    let marks: number = MARK_NONE;
    if (chars === 0 && linear && v.height >= 4) {
      for (let e = win[0]; e < win[1]; e++) {
        if (o.main[e] === 0) continue;
        if (ofColor(v, o, b, e, s.map) < 0) {
          mk[0] = room[0];
          marks = markForm(v, y0 | 0, y1 | 0);
        }
        break;
      }
    }
    if (thin && chars === 0 && marks === MARK_NONE) continue;
    let edge = 0;
    if (chars > 0) {
      letterFont(v, font);
      letterLines(v, font, y0, y1, baselines);
      monoCell(tm, font);
      edge = (chars * cellW[0] + (chars === 3 ? 2 * OF_MARGIN : 0)) * dpr;
    }
    for (let e = win[0]; e < win[1]; e++) {
      if (o.stop[e] <= hull[0]) continue;
      const n = placeRun(v, b, o.start[e], o.stop[e], piece);
      if (n === 0) continue;
      const color = ofColor(v, o, b, e, s.map);
      if (!thin && color >= 0) {
        rects.begin(color, ofCells.alpha(b, o.start[e]));
        for (let q = 0; q < n; q++) rects.push(piece[2 * q], y0, Math.max(1, piece[2 * q + 1] - piece[2 * q]), y1 - y0);
      }
      if (o.main[e] === 0 || (chars === 0 && (marks === MARK_NONE || color >= 0))) continue;
      let q = 0;
      for (let r = 1; r < n; r++) if (piece[2 * r + 1] - piece[2 * r] > piece[2 * q + 1] - piece[2 * q]) q = r;
      const x0 = piece[2 * q];
      const x1 = piece[2 * q + 1];
      if (chars === 0) {
        if (marks === MARK_DOTS) pushDot(x0 | 0, x1 | 0);
        else pushBar(x0 | 0, x1 | 0);
        continue;
      }
      if (linear && span > 1 && x1 - x0 >= 4) {
        lines.begin(SkinToken.surface);
        lines.push(x0, y0, 1, y1 - y0);
      }
      // A residue cut by the view's edge has a narrower part than the row's text, with a three-letter code's margins
      // (a device px of rounding allowed): its text would reach into its neighbor, so none.
      if (x1 - x0 + 1 < edge) continue;
      // A residue the map changes is drawn as the one letter shown (its three-letter form is the stored one's).
      const c0 = o.code[e];
      const c = s.map !== null ? s.map[c0 & 127] : c0;
      letters.begin(color >= 0 ? cellInk(color) : SkinToken.ink);
      const long = chars === 3 && c === c0;
      const t = letters.add(long ? three(c) : CHARS[c & 127]);
      letters.x[t] = (x0 + x1) / 2;
      letters.y[t] = baselines[long ? 0 : LETTER_LINE[c & 127]];
    }
    if (marks === MARK_BARS) closeBar();
  }
}

/** Below the row's minimum: the same draw (cells only: letters need a row to read them in). */
export const LettersRenderer = kernelEntry<LettersState>('letters', 22, {
  geometries: ['linear', 'wrapped'],
  stable,
  miniature: false,
  hit,
  draw,
  tiny: draw,
});
