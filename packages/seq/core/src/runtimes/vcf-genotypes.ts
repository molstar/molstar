/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { categoryColumn } from '@molstar/core/repr/tree';
import {
  type MatrixItems,
  type MatrixSummary,
  type RegionAnswer,
  type RegionRequest,
  VcfFixedFieldNames,
} from '../data/protocol.js';
import { failure } from './internal/answers.js';
import { floatOf, type VcfRows } from './vcf-parse.js';
import type { GenotypeSchema, MatrixBounds } from './vcf-samples.js';

// A VCF's matrix (_spec/README.md "Data", VCF: the matrix): items the ALT alleles of a region × states the samples
// held (the header's order, narrowed by the source's `samples`), `codes` the allele's genotype class in `GT` (0 no copy,
// 1 some but not every copy, 2 every copy (a haploid `1` too), 3 partly called, 255 uncalled), or `value` float32 by
// `RegionRequest.field` (`GT` the allele's copies; a numeric FORMAT key of `Number=1` or `A`; NaN missing). Genotypes
// are read only for a matrix request and only for the samples held, from the lines' bytes, once a line; a coarser
// level is a summary counted from the region's items (majority class and its share, or mean, min, max and count).
// Loaded on the first matrix request (`vcf-answer.ts`), so a host reading sites never loads it.

const tab = 9;
const colon = 58;
const comma = 44;
const dot = 46;

/** The fixed fields' places among a row's values: the record's ID, the allele's category. */
const Label = VcfFixedFieldNames.indexOf('label');
const Category = VcfFixedFieldNames.indexOf('category');

/** The longest call kept as text (`0|1|1|0`): a longer one is read as its class alone. */
const MaxCallBytes = 16;

/**
 * The distinct genotype calls of one matrix read, by their bytes (a hash, checked against the text): an index a cell
 * holds instead of a string. At most 255 distinct; a call past that, or past `MaxCallBytes`, has index 0 (its class).
 */
class CallDictionary {
  readonly texts: string[] = [''];
  private readonly byHash = new Map<number, number>();

  index(bytes: Uint8Array, s: number, e: number): number {
    const n = e - s;
    if (n <= 0 || n > MaxCallBytes) return 0;
    let h = n;
    for (let i = s; i < e; i++) h = Math.imul(h ^ bytes[i], 16777619);
    const k = this.byHash.get(h);
    if (k !== undefined) {
      const t = this.texts[k];
      if (t.length !== n) return 0;
      for (let i = 0; i < n; i++) if (t.charCodeAt(i) !== bytes[s + i]) return 0;
      return k;
    }
    if (this.texts.length >= 256) return 0;
    let text = '';
    for (let i = s; i < e; i++) text += String.fromCharCode(bytes[i]);
    this.texts.push(text);
    this.byHash.set(h, this.texts.length - 1);
    return this.texts.length - 1;
  }
}

/**
 * Buffers one matrix read reuses across its lines, sized once for its states and its highest allele (counted in its
 * bound): per state the ploidy, the `.` alleles and the field's span, per allele and state the copies.
 */
class LineGenotypes {
  /** Per state, the call's index in the read's `CallDictionary` (0 none). */
  readonly call: Uint8Array;
  readonly ploidy: Uint8Array;
  readonly dots: Uint8Array;
  readonly copies: Uint8Array;
  readonly fieldStart: Int32Array;
  readonly fieldEnd: Int32Array;

  constructor(
    states: number,
    alleles: number,
    readonly calls = new CallDictionary(),
  ) {
    this.call = new Uint8Array(states);
    this.ploidy = new Uint8Array(states);
    this.dots = new Uint8Array(states);
    this.fieldStart = new Int32Array(states);
    this.fieldEnd = new Int32Array(states);
    this.copies = new Uint8Array(states * alleles);
  }

  reset(states: number, alleles: number) {
    this.call.fill(0, 0, states);
    this.ploidy.fill(0, 0, states);
    this.dots.fill(0, 0, states);
    this.fieldStart.fill(-1, 0, states);
    this.fieldEnd.fill(-1, 0, states);
    this.copies.fill(0, 0, states * alleles);
  }
}

/** The position of the `n`th tab of `bytes[s, e)` plus one (the column after it), or −1. */
function afterTabs(bytes: Uint8Array, s: number, e: number, n: number): number {
  let at = s;
  for (let k = 0; k < n; k++) {
    const t = bytes.indexOf(tab, at);
    if (t < 0 || t >= e) return -1;
    at = t + 1;
  }
  return at;
}

/** The index of `key` among a FORMAT column's `:`-separated keys `bytes[s, e)`, or −1. */
function formatIndex(bytes: Uint8Array, s: number, e: number, key: string): number {
  let k = 0;
  for (let a = s; a <= e; k++) {
    let b = a;
    while (b < e && bytes[b] !== colon) b++;
    if (b - a === key.length) {
      let same = true;
      for (let i = 0; i < key.length && same; i++) same = bytes[a + i] === key.charCodeAt(i);
      if (same) return k;
    }
    a = b + 1;
  }
  return -1;
}

/** The span of a sample column's `index`th `:` subfield in `bytes[s, e)`, written to `out`; false when it has none. */
function subfield(bytes: Uint8Array, s: number, e: number, index: number, out: Int32Array): boolean {
  let a = s;
  for (let k = 0; k < index; k++) {
    while (a < e && bytes[a] !== colon) a++;
    if (a >= e) return false;
    a++;
  }
  let b = a;
  while (b < e && bytes[b] !== colon) b++;
  out[0] = a;
  out[1] = b;
  return true;
}

const span = new Int32Array(2);

/**
 * Reads a line's genotypes for the held states (`columns`): each state's ploidy, `.` alleles and the copies of each ALT
 * (`gt` the GT key's index), and the span of the value field (`fi`, its index; −1 none). Columns past the line are
 * uncalled. One pass over the line's sample columns.
 */
function readLine(
  bytes: Uint8Array,
  s: number,
  e: number,
  columns: Int32Array,
  alleles: number,
  g: LineGenotypes,
  field: string | undefined,
) {
  const S = columns.length;
  g.reset(S, alleles);
  const fmt = afterTabs(bytes, s, e, 8);
  if (fmt < 0) return;
  let fmtEnd = bytes.indexOf(tab, fmt);
  if (fmtEnd < 0 || fmtEnd > e) fmtEnd = e;
  const gt = formatIndex(bytes, fmt, fmtEnd, 'GT');
  const fi = field === undefined || field === 'GT' ? -1 : formatIndex(bytes, fmt, fmtEnd, field);
  let at = fmtEnd + 1;
  let col = 0;
  for (let st = 0; st < S && at <= e; st++) {
    const want = columns[st];
    // the columns before this state's read past, a tab each
    while (col < want && at <= e) {
      const t = bytes.indexOf(tab, at);
      at = t < 0 || t >= e ? e + 1 : t + 1;
      col++;
    }
    if (at > e) break;
    let end = bytes.indexOf(tab, at);
    if (end < 0 || end > e) end = e;
    if (gt >= 0 && subfield(bytes, at, end, gt, span)) {
      let p = 0;
      let d = 0;
      for (let i = span[0]; i < span[1]; ) {
        const c = bytes[i];
        if (c === dot) {
          d++;
          p++;
          i++;
        } else if (c >= 48 && c <= 57) {
          let v = 0;
          while (i < span[1] && bytes[i] >= 48 && bytes[i] <= 57 && v < 65536) v = v * 10 + bytes[i++] - 48;
          p++;
          if (v > 0 && v <= alleles) {
            const x = (v - 1) * S + st;
            if (g.copies[x] < 255) g.copies[x]++;
          }
        } else i++;
      }
      g.call[st] = g.calls.index(bytes, span[0], span[1]);
      g.ploidy[st] = Math.min(255, p);
      g.dots[st] = Math.min(255, d);
    }
    if (fi >= 0 && subfield(bytes, at, end, fi, span)) {
      g.fieldStart[st] = span[0];
      g.fieldEnd[st] = span[1];
    }
    at = end + 1;
    col++;
  }
}

/** The genotype class of an allele's copies in a call of `ploidy` alleles, `dots` of them missing. */
export function genotypeClass(copies: number, ploidy: number, dots: number): number {
  if (ploidy === 0 || dots === ploidy) return 255;
  if (dots > 0) return 3;
  if (copies === 0) return 0;
  return copies === ploidy ? 2 : 1;
}

/** A decimal number in `bytes[s, e)` (`12`, `-0.5`, `1e-3`), NaN missing (`.`, empty) or not a number. */
function numberAt(bytes: Uint8Array, s: number, e: number): number {
  if (e <= s || (e - s === 1 && bytes[s] === dot)) return Number.NaN;
  let i = s;
  let sign = 1;
  if (bytes[i] === 45 || bytes[i] === 43) sign = bytes[i++] === 45 ? -1 : 1;
  // digits, a fraction and an exponent read from the bytes (15 significant digits kept, float32 holds fewer)
  let v = 0;
  let digits = 0;
  let scale = 0;
  while (i < e && bytes[i] >= 48 && bytes[i] <= 57) {
    if (digits < 15) v = v * 10 + bytes[i] - 48;
    else scale++;
    digits++;
    i++;
  }
  if (i < e && bytes[i] === dot) {
    i++;
    while (i < e && bytes[i] >= 48 && bytes[i] <= 57) {
      if (digits < 15) {
        v = v * 10 + bytes[i] - 48;
        scale--;
      }
      digits++;
      i++;
    }
  }
  if (!digits) return wordAt(bytes, s, e);
  if (i < e && (bytes[i] === 101 || bytes[i] === 69)) {
    i++;
    let es = 1;
    if (bytes[i] === 45 || bytes[i] === 43) es = bytes[i++] === 45 ? -1 : 1;
    let x = 0;
    const from = i;
    while (i < e && bytes[i] >= 48 && bytes[i] <= 57 && x < 1000) x = x * 10 + bytes[i++] - 48;
    if (i === from) return Number.NaN;
    scale += es * x;
  }
  if (i !== e) return Number.NaN;
  return sign * (scale < 0 ? v / 10 ** -scale : v * 10 ** scale);
}

/** `Inf`, `Infinity` or `NaN` (any case, a sign): VCF's words for a float; else NaN. Only a cell with no digit gets here. */
function wordAt(bytes: Uint8Array, s: number, e: number): number {
  if (e - s > 9) return Number.NaN;
  let text = '';
  for (let k = s; k < e; k++) text += String.fromCharCode(bytes[k]);
  return floatOf(text);
}

/** The `k`th (0-based) `,`-separated part of `bytes[s, e)` as a number, NaN when absent. */
function partAt(bytes: Uint8Array, s: number, e: number, k: number): number {
  let a = s;
  for (let i = 0; i < k; i++) {
    while (a < e && bytes[a] !== comma) a++;
    if (a >= e) return Number.NaN;
    a++;
  }
  let b = a;
  while (b < e && bytes[b] !== comma) b++;
  return numberAt(bytes, a, b);
}

/**
 * The matrix of `rows` (row numbers in region order) over `request`: at level 1, or as items (`request.items`), the
 * items × states; coarser, a summary counted from them. Past `maxBytes` (counted before it is allocated) `too-large`.
 */
export function matrixAnswer(
  data: VcfRows | undefined,
  rows: readonly number[],
  g: GenotypeSchema,
  request: RegionRequest,
  bounds: MatrixBounds,
): RegionAnswer {
  const S = g.states.length;
  const field = request.field;
  if (field !== undefined && !g.number.has(field)) return failure('not-served', `No genotype field ${field}`);
  const values = field !== undefined;
  const summary = request.level > 1 && !request.items;
  const n = rows.length;
  const L = request.level;
  const lo = Math.floor(request.start / L);
  const bins = summary ? Math.max(0, Math.ceil(request.end / L) - lo) : 0;
  // the highest allele read: a line's copies are counted a state and an allele
  let alleles = 0;
  if (data) for (let k = 0; k < n; k++) alleles = Math.max(alleles, data.allele[rows[k]]);
  // the items' cells, a summary's with what it counts in (codes: four counts and two bytes; values: five and a sum),
  // and the line's buffers (copies a state and an allele; ploidy, missing alleles and a field's span a state)
  const bytes = n * S * (values ? 4 : 2) + (summary ? bins * S * (values ? 24 : 18) : 0) + (n ? S * (alleles + 10) : 0);
  if (bytes > bounds.maxBytes)
    return failure('too-large', `A matrix of ${n} items × ${S} samples: ${bytes} bytes, over ${bounds.maxBytes}`);
  const genotypes = new LineGenotypes(S, alleles);
  const codes = values ? undefined : new Uint8Array(n * S);
  const calls = values ? undefined : new Uint8Array(n * S);
  const value = values ? new Float32Array(n * S) : undefined;
  const number = values ? g.number.get(field!)! : '1';
  if (data) {
    let line = -1;
    for (let i = 0; i < n; i++) {
      const r = rows[i];
      const s = data.lineStart[r];
      // rows of one line read it once (its alleles are its rows)
      if (s !== line) {
        line = s;
        // a line's rows are next to each other (in order of start, then of the file): its alleles, read forward
        let alleles = 0;
        for (let k = i; k < n && data.lineStart[rows[k]] === s; k++) alleles = Math.max(alleles, data.allele[rows[k]]);
        readLine(data.bytes, s, data.lineEnd[r], g.columns, alleles, genotypes, field);
      }
      const a = data.allele[r] - 1;
      for (let st = 0; st < S; st++) {
        const c = genotypes.copies[a * S + st];
        const p = genotypes.ploidy[st];
        const d = genotypes.dots[st];
        if (codes) {
          codes[i * S + st] = genotypeClass(c, p, d);
          calls![i * S + st] = genotypes.call[st];
        } else if (field === 'GT') value![i * S + st] = p === 0 || d === p ? Number.NaN : c;
        else {
          const fs = genotypes.fieldStart[st];
          value![i * S + st] =
            fs < 0
              ? Number.NaN
              : number === 'A'
                ? partAt(data.bytes, fs, genotypes.fieldEnd[st], a)
                : numberAt(data.bytes, fs, genotypes.fieldEnd[st]);
        }
      }
    }
  }
  const start = data ? Int32Array.from(rows, (r) => data.start[r]) : new Int32Array(0);
  const end = data ? Int32Array.from(rows, (r) => data.end[r]) : new Int32Array(0);
  if (!summary) {
    const items: MatrixItems = {
      form: 'items',
      start,
      end,
      states: S,
      ...(codes ? { codes, calls: calls!, callTexts: genotypes.calls.texts } : { value: value! }),
      id: categoryColumn(data ? rows.map((r) => data.key[r]) : []),
      // the record's ID and the allele's category (the fixed fields' first and last): the cell's title
      label: categoryColumn(data ? rows.map((r) => data.values[Label][r] as string | null) : []),
      category: categoryColumn(data ? rows.map((r) => data.values[Category][r] as string | null) : []),
    };
    return { ok: true, tables: { matrix: items } };
  }
  return { ok: true, tables: { matrix: summarize(start, end, codes, value, S, lo, bins, L, request) } };
}

/**
 * A summary of items × states at `level` over the bins `[lo, lo + bins)`: an item counted in every bin its extent
 * touches inside the region (as a ranges density); codes as each bin's majority class (0 to 3, ties to the lower) and
 * its share of the called (0 to 255; 255 and 0 where none is called), values as mean, min, max and count (NaN not
 * counted).
 */
function summarize(
  start: Int32Array,
  end: Int32Array,
  codes: Uint8Array | undefined,
  value: Float32Array | undefined,
  S: number,
  lo: number,
  bins: number,
  level: number,
  request: RegionRequest,
): MatrixSummary {
  const cells = bins * S;
  const binsOf = (k: number): [number, number] => {
    const a = Math.max(lo, Math.floor(Math.max(start[k], request.start) / level));
    const b = Math.min(lo + bins - 1, Math.floor((Math.min(Math.max(end[k], start[k] + 1), request.end) - 1) / level));
    return [a - lo, b - lo];
  };
  if (codes) {
    const counts = new Float32Array(cells * 4);
    for (let k = 0; k < start.length; k++) {
      const [a, b] = binsOf(k);
      for (let x = a; x <= b; x++)
        for (let s = 0; s < S; s++) {
          const c = codes[k * S + s];
          if (c <= 3) counts[(x * S + s) * 4 + c]++;
        }
    }
    const classes = new Uint8Array(cells).fill(255);
    const shares = new Uint8Array(cells);
    for (let i = 0; i < cells; i++) {
      const called = counts[i * 4] + counts[i * 4 + 1] + counts[i * 4 + 2] + counts[i * 4 + 3];
      if (!called) continue;
      let best = 0;
      for (let c = 1; c < 4; c++) if (counts[i * 4 + c] > counts[i * 4 + best]) best = c;
      classes[i] = best;
      shares[i] = Math.round((255 * counts[i * 4 + best]) / called);
    }
    return { form: 'summary', start: lo * level, states: S, classes, shares };
  }
  const mean = new Float32Array(cells);
  const min = new Float32Array(cells).fill(Number.NaN);
  const max = new Float32Array(cells).fill(Number.NaN);
  const count = new Float32Array(cells);
  const sum = new Float64Array(cells);
  for (let k = 0; k < start.length; k++) {
    const [a, b] = binsOf(k);
    for (let x = a; x <= b; x++)
      for (let s = 0; s < S; s++) {
        const v = value![k * S + s];
        if (Number.isNaN(v)) continue;
        const i = x * S + s;
        sum[i] += v;
        count[i]++;
        if (!(v >= min[i])) min[i] = v;
        if (!(v <= max[i])) max[i] = v;
      }
  }
  for (let i = 0; i < cells; i++) mean[i] = count[i] ? sum[i] / count[i] : Number.NaN;
  return { form: 'summary', start: lo * level, states: S, mean, min, max, count };
}
