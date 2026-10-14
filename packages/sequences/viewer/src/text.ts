/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Text measuring (seqspec:viewer#text-measuring): prepare measures once per font into a cache; layout is arithmetic
// over the
// prepared arrays, with no measuring, no strings and no allocation. Labels, tick numbers and letter cells here;
// wrapped runs in runs.ts.
//
// Adapted from pretext (https://github.com/chenglou/pretext, commit e699e27e7e4b48a50024fca545c6a12bc6ac975d):
// the prepare/layout split, a width cache per font dropped when fonts change (`clearCache`), per-grapheme advances
// measured as prefixes to cut inside a segment (its `segment-prefixes` fit), graphemes only where text is not ASCII.
// Not taken: its whitespace modes, line-breaking rules per script and engine, bidi and emoji correction
// (README.md open question 3).
//
// pretext is MIT licensed: Copyright (c) 2026 Pretext contributors. Permission is hereby granted, free of charge,
// to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in
// the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions: The above copyright notice and this permission notice
// shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import type { Alphabet } from '@molstar/seq-core/model';
import { type SeqState, state } from './state.js';
import { hasMetrics } from './svg/metrics.js';

/** `ports.text`: a width in CSS px. The engine's is canvas `measureText`. */
export interface TextMeasurer {
  measure(font: string, text: string): number;
}
/** An interned font string; one cache per font. */
export type FontId = number;

export interface PreparedLabel {
  readonly text: string;
  readonly width: number;
  /** Cumulative advance at the end of each grapheme. */
  readonly ends: Float32Array;
  /** Each grapheme's start in the string, and the string's length last; null: ASCII, one unit each. */
  readonly offsets: Uint16Array | null;
}

/** The px of a CSS font string (`"600 12px sans-serif"`); 0 when none is given. */
export function fontPx(font: string): number {
  const m = /([0-9.]+)px/.exec(font);
  return m ? Number(m[1]) : 0;
}

/** The default port (seqspec:viewer#the-context): a monospace estimate, 0.6 em per code point. */
export const monospaceMeasurer: TextMeasurer = {
  measure: (font, text) => {
    let n = 0;
    for (let i = 0; i < text.length; i++) if ((text.charCodeAt(i) & 0xfc00) !== 0xdc00) n++;
    return 0.6 * fontPx(font) * n;
  },
};

/** A capital's and a lower-case letter's height in em where the measurer does not know them (`ink`). */
const CAP_EM = 0.72;
const X_EM = 0.54;

const LETTERS: Record<Alphabet, string> = {
  protein: 'ACDEFGHIKLMNPQRSTVWYBZJUOX',
  dna: 'ACGTN',
  rna: 'ACGUN',
  generic: 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',
};

class Label implements PreparedLabel {
  readonly text: string;
  readonly width: number;
  readonly ends: Float32Array;
  readonly offsets: Uint16Array | null;
  /** Its insertion stamp in the font's cache: a hit on an old one moves it to the end (an LRU, edge-rate). */
  stamp: number;
  /** The strings drawn, per grapheme count. */
  cut: Map<number, string> | null = null;

  constructor(text: string, ends: Float32Array, offsets: Uint16Array | null, stamp: number) {
    this.text = text;
    this.ends = ends;
    this.offsets = offsets;
    this.width = ends.length ? ends[ends.length - 1] : 0;
    this.stamp = stamp;
  }
}

interface FontCache {
  font: string;
  readonly role: 'mono' | 'ui';
  readonly px: number;
  readonly weight: number;
  advances: Float32Array | null;
  readonly cells: Map<string, number>;
  readonly labels: Map<string, Label>;
  /** `ink` by DPR × 64. */
  readonly ink: Map<number, Float64Array>;
  stamp: number;
}

export interface TextManagerApi {
  /** The measurer in use (`ports.text`, or the one `setTextMeasurer` gave). */
  readonly port: TextMeasurer;
  /** A new measurer or a loaded font: every prepared text is stale. */
  readonly version: SeqState<number>;
  /** Interned by (role, px to a quarter, weight): the table stays bounded whatever sizes are asked for. */
  font(role: 'mono' | 'ui', px: number, weight?: 400 | 600): FontId;
  /** The CSS font string a canvas draws it with. */
  fontString(font: FontId): string;
  advances(font: FontId): Float32Array;
  /**
   * The ink heights of a capital ("H") and of a lower-case letter ("x") above the alphabetic baseline, in device px at
   * a DPR (`dpr64`: DPR × 64, an integer, so nothing boxes): from the measurer's `metrics` at the device size (the
   * canvas measurer's are read from its raster, hinting included), else 0.72 and 0.54 em. Letters are centered on them
   * (text placement). Read once per font and DPR.
   */
  ink(font: FontId, dpr64: number): Float64Array;
  /** The widest letter of an alphabet: the letters LOD test. */
  cell(font: FontId, alphabet: Alphabet | null): number;
  /** Digits, sign and suffix from advances, in whole px rounded up; no string. */
  numberWidth(font: FontId, n: number, suffix?: number): number;
  /** Measured once per (text, font); an LRU of `Text.CacheEntries` per font. */
  label(text: string, font: FontId): PreparedLabel;
  /** The families fonts are made of (`setFamilies`). */
  readonly families: { readonly mono: string; readonly ui: string };
  /** The skin's families (the engine's skin reader): a change remakes every font string and drops every cache. */
  setFamilies(mono: string, ui: string): void;
  /** A web font loaded after text was prepared: every cache dropped, `version` bumped; the measurer stays. */
  invalidate(): void;
}

/** `ctx.text` (seqspec:viewer#state): fonts interned, advances and labels cached per font, all dropped together. */
export class TextManager implements TextManagerApi {
  readonly version = state(0);
  /** The families fonts are made of; the engine sets the skin's (R1). */
  readonly families = { mono: 'monospace', ui: 'sans-serif' };
  private measurer: TextMeasurer;
  private readonly capacity: () => number;
  private readonly ids = new Map<number, FontId>();
  private fonts: FontCache[] = [];
  private segmenter: Intl.Segmenter | null = null;

  constructor(measurer: TextMeasurer, capacity: () => number) {
    this.measurer = measurer;
    this.capacity = capacity;
  }

  get port(): TextMeasurer {
    return this.measurer;
  }
  setMeasurer(m: TextMeasurer) {
    this.measurer = m;
    this.invalidate();
  }
  setFamilies(mono: string, ui: string) {
    if (mono === this.families.mono && ui === this.families.ui) return;
    this.families.mono = mono;
    this.families.ui = ui;
    for (const f of this.fonts) f.font = `${f.weight} ${f.px}px ${this.families[f.role]}`;
    this.invalidate();
  }
  /** A web font loaded after text was prepared (pretext's `clearCache`): every cache is dropped. */
  invalidate() {
    this.fonts = this.fonts.map((f) => ({
      ...f,
      advances: null,
      cells: new Map(),
      labels: new Map(),
      ink: new Map(),
      stamp: 0,
    }));
    this.version.set(this.version.value + 1);
  }

  font(role: 'mono' | 'ui', px0: number, weight: 400 | 600 = 400): FontId {
    // A quarter px at most: a size that follows a zoom continuously makes a few fonts a px, never one a frame (each
    // a font string, its advances measured, its ink read and its atlas cells).
    const px = Math.round(px0 * 4) / 4;
    const key = (px * 1024 + weight) * 2 + (role === 'mono' ? 1 : 0);
    const id = this.ids.get(key);
    if (id !== undefined) return id;
    const font = `${weight} ${px}px ${this.families[role]}`;
    this.fonts.push({
      font,
      role,
      px,
      weight,
      advances: null,
      cells: new Map(),
      labels: new Map(),
      ink: new Map(),
      stamp: 0,
    });
    this.ids.set(key, this.fonts.length - 1);
    return this.fonts.length - 1;
  }
  fontString(font: FontId): string {
    return this.fonts[font].font;
  }

  advances(font: FontId): Float32Array {
    const f = this.fonts[font];
    if (f.advances) return f.advances;
    const a = new Float32Array(128);
    for (let c = 32; c < 127; c++) a[c] = this.measurer.measure(f.font, String.fromCharCode(c));
    f.advances = a;
    return a;
  }

  ink(font: FontId, dpr64: number): Float64Array {
    const f = this.fonts[font];
    let a = f.ink.get(dpr64);
    if (a) return a;
    const px = (f.px * dpr64) / 64;
    const m = this.measurer;
    const box = hasMetrics(m) ? m.metrics(`${f.weight} ${px}px ${this.families[f.role]}`) : null;
    a = new Float64Array(2);
    a[0] = (box?.cap ?? CAP_EM) * px;
    a[1] = (box?.xHeight ?? X_EM) * px;
    f.ink.set(dpr64, a);
    return a;
  }

  cell(font: FontId, alphabet: Alphabet | null): number {
    const f = this.fonts[font];
    const key = alphabet ?? 'generic';
    let w = f.cells.get(key);
    if (w === undefined) {
      const a = this.advances(font);
      w = 0;
      for (const c of LETTERS[key]) w = Math.max(w, a[c.charCodeAt(0)]);
      f.cells.set(key, w);
    }
    return w;
  }

  numberWidth(font: FontId, n: number, suffix = -1): number {
    const a = this.fonts[font].advances ?? this.advances(font);
    let v = n < 0 ? -n : n;
    let w = n < 0 ? a[45] : 0;
    do {
      const q = (v / 10) | 0;
      w += a[48 + v - q * 10];
      v = q;
    } while (v > 0);
    // Whole px, rounded up (enough to thin ticks): a small integer, so the result never boxes at the call.
    return Math.ceil(suffix >= 0 && suffix < 128 ? w + a[suffix] : w);
  }

  label(text: string, font: FontId): PreparedLabel {
    const f = this.fonts[font];
    const hit = f.labels.get(text);
    const cap = Math.max(1, this.capacity());
    if (hit) {
      if (f.stamp - hit.stamp > cap / 2) {
        f.labels.delete(text);
        hit.stamp = ++f.stamp;
        f.labels.set(text, hit);
      }
      return hit;
    }
    const label = this.prepare(text, f, font);
    f.labels.set(text, label);
    if (f.labels.size > cap) f.labels.delete(f.labels.keys().next().value as string);
    return label;
  }

  private prepare(text: string, f: FontCache, font: FontId): Label {
    let ascii = true;
    for (let i = 0; i < text.length && ascii; i++) ascii = text.charCodeAt(i) < 128;
    if (ascii) {
      const a = this.advances(font);
      const ends = new Float32Array(text.length);
      let w = 0;
      for (let i = 0; i < text.length; i++) ends[i] = w += a[text.charCodeAt(i)];
      return new Label(text, ends, null, ++f.stamp);
    }
    // Graphemes (Intl.Segmenter), each cut measured as the prefix that ends there, as pretext fits a segment.
    this.segmenter ??= new Intl.Segmenter(undefined, { granularity: 'grapheme' });
    const starts: number[] = [];
    for (const g of this.segmenter.segment(text)) starts.push(g.index);
    const offsets = Uint16Array.from([...starts, text.length]);
    const ends = new Float32Array(starts.length);
    for (let g = 0; g < starts.length; g++) ends[g] = this.measurer.measure(f.font, text.slice(0, offsets[g + 1]));
    return new Label(text, ends, offsets, ++f.stamp);
  }
}

/** The graphemes that fit `maxWidth`, leaving room for an ellipsis when cut; a binary search over `ends`. */
export function fitCount(p: PreparedLabel, maxWidth: number, ellipsis: number): number {
  fitArgs[0] = maxWidth;
  fitArgs[1] = ellipsis;
  return fitIn(p, fitArgs);
}
const fitArgs = new Float64Array(2);
/** `fitCount` with (maxWidth, ellipsis) in `args`: per-frame callers pass no double through a call (V8 boxes it). */
export function fitIn(p: PreparedLabel, args: Float64Array): number {
  const maxWidth = args[0];
  const ellipsis = args[1];
  const ends = p.ends;
  const n = ends.length;
  if (n === 0 || p.width <= maxWidth) return n;
  let lo = 0;
  let hi = n - 1;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (ends[mid - 1] + ellipsis <= maxWidth) lo = mid;
    else hi = mid - 1;
  }
  return lo;
}

/** The string to draw for `count` graphemes (an ellipsis when cut), built once per (label, count). */
export function fitted(p: PreparedLabel, count: number): string {
  const n = p.ends.length;
  if (count >= n) return p.text;
  const l = p as Label;
  l.cut ??= new Map();
  let s = l.cut.get(count);
  if (s === undefined) {
    s = `${p.text.slice(0, p.offsets ? p.offsets[count] : count)}…`;
    l.cut.set(count, s);
  }
  return s;
}
