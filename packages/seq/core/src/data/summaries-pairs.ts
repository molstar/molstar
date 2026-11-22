/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, categoryFromCodes, type NumericColumn } from '@molstar/core/repr/tree';
import type { PairsBins, PairsList } from './protocol.js';

// The pairs level (_spec/README.md "Data", Summaries): bin-pair counts at a level, by relative strand and by `type`
// when the coverage lists it (`binBy`), a field summed when asked. `binPairs` makes a level from items, `coarsenPairBins`
// a coarser level from bins (exact: levels are nested, bins aligned on position 0). An end's bin is its first position's.

/** What `binPairs` keeps besides the bin pair: the `type` (`category`), a numeric field summed (`field`). */
export interface BinPairsOptions {
  readonly category?: boolean;
  readonly field?: string;
  /** Both ends' system when the list carries none (written into the bins when `systems` is asked). */
  readonly system?: string;
  /** Write each end's system even when the list has none (a pyramid over several systems). */
  readonly systems?: boolean;
}

/** Category codes of `c` into one dictionary shared with others (`dict`, `index`): −1 stays missing. */
function remap(c: CategoryColumn | undefined, dict: string[], index: Map<string, number>): Int32Array | undefined {
  if (!c) return undefined;
  const local = c.dict.map((v) => {
    let i = index.get(v);
    if (i === undefined) {
      i = dict.length;
      dict.push(v);
      index.set(v, i);
    }
    return i;
  });
  const out = new Int32Array(c.codes.length);
  for (let i = 0; i < out.length; i++) out[i] = c.codes[i] < 0 ? -1 : local[c.codes[i]];
  return out;
}

/** Grows a typed array to hold `n`, doubling. */
function grow<A extends Int32Array | Float64Array>(a: A, n: number): A {
  if (n <= a.length) return a;
  const out = new (a.constructor as new (n: number) => A)(Math.max(n, 2 * a.length));
  out.set(a);
  return out;
}

/**
 * Bins under construction: each bin pair's group (both systems, strand, category), bins, count and sum, found by an
 * open-addressing table over typed arrays (a pyramid of 600,000 pairs is a few dozen slices: no map a pair).
 */
export class Binner {
  n = 0;
  aSys = new Int32Array(64);
  bSys = new Int32Array(64);
  aBin = new Int32Array(64);
  bBin = new Int32Array(64);
  strand = new Int32Array(64);
  cat = new Int32Array(64);
  count = new Float64Array(64);
  sum = new Float64Array(64);
  /** The groups seen (small: systems² × strands × categories), each an id. */
  private readonly _groups = new Map<number | string, number>();
  private readonly _group: number[] = [];
  private _lastKey: number | string | undefined;
  private _lastGroup = -1;
  /** Slots: a row + 1, 0 empty. */
  private _slots = new Int32Array(1024);

  private _groupOf(as: number, bs: number, st: number, cat: number): number {
    const key =
      as < 0xffff && bs < 0xffff && cat < 0xffff
        ? (((as + 1) * 0x10000 + (bs + 1)) * 4 + (st + 1)) * 0x10001 + (cat + 1)
        : `${as},${bs},${st},${cat}`;
    if (key === this._lastKey) return this._lastGroup;
    let g = this._groups.get(key);
    if (g === undefined) this._groups.set(key, (g = this._groups.size));
    this._lastKey = key;
    this._lastGroup = g;
    return g;
  }

  private _hash(g: number, ab: number, bb: number): number {
    let h = Math.imul(ab, 0x9e3779b1) ^ Math.imul(bb + 0x7f4a7c15, 0x85ebca6b) ^ Math.imul(g + 1, 0xc2b2ae35);
    h ^= h >>> 15;
    return h & (this._slots.length - 1);
  }

  add(as: number, bs: number, ab: number, bb: number, st: number, cat: number, n: number, s: number): void {
    const g = this._groupOf(as, bs, st, cat);
    const slots = this._slots;
    const mask = slots.length - 1;
    let h = this._hash(g, ab, bb);
    for (;;) {
      const at = slots[h];
      if (at === 0) break;
      const r = at - 1;
      if (this.aBin[r] === ab && this.bBin[r] === bb && this._group[r] === g) {
        this.count[r] += n;
        this.sum[r] += s;
        return;
      }
      h = (h + 1) & mask;
    }
    const r = this.n++;
    if (r >= this.aBin.length) {
      const m = r + 1;
      this.aSys = grow(this.aSys, m);
      this.bSys = grow(this.bSys, m);
      this.aBin = grow(this.aBin, m);
      this.bBin = grow(this.bBin, m);
      this.strand = grow(this.strand, m);
      this.cat = grow(this.cat, m);
      this.count = grow(this.count, m);
      this.sum = grow(this.sum, m);
    }
    this.aSys[r] = as;
    this.bSys[r] = bs;
    this.aBin[r] = ab;
    this.bBin[r] = bb;
    this.strand[r] = st;
    this.cat[r] = cat;
    this.count[r] = n;
    this.sum[r] = s;
    this._group[r] = g;
    slots[h] = r + 1;
    // at half full, twice the slots
    if (2 * this.n > slots.length) this._rehash();
  }

  private _rehash(): void {
    const slots = new Int32Array(this._slots.length * 2);
    this._slots = slots;
    const mask = slots.length - 1;
    for (let r = 0; r < this.n; r++) {
      let h = this._hash(this._group[r], this.aBin[r], this.bBin[r]);
      while (slots[h] !== 0) h = (h + 1) & mask;
      slots[h] = r + 1;
    }
  }

  finish(
    sysDict: readonly string[] | undefined,
    catDict: readonly string[] | undefined,
    strand: boolean,
    sum: boolean,
  ) {
    const n = this.n;
    const out: { -readonly [K in keyof PairsBins]: PairsBins[K] } = {
      form: 'bins',
      aBin: this.aBin.slice(0, n),
      bBin: this.bBin.slice(0, n),
      count: Float32Array.from(this.count.subarray(0, n)),
    };
    if (sysDict) {
      out.aSystem = categoryFromCodes(sysDict, this.aSys.subarray(0, n));
      out.bSystem = categoryFromCodes(sysDict, this.bSys.subarray(0, n));
    }
    if (strand) out.strand = Int8Array.from(this.strand.subarray(0, n));
    if (catDict) out.category = categoryFromCodes(catDict, this.cat.subarray(0, n));
    if (sum) out.sum = Float32Array.from(this.sum.subarray(0, n));
    return out as PairsBins;
  }
}

/** A field's numeric column of a list, when it has one by that name. */
function numericField(t: PairsList, name: string | undefined): NumericColumn | undefined {
  const c = name === undefined ? undefined : t.fields?.find((f) => f.name === name)?.column;
  return c && !(c instanceof CategoryColumn) ? c : undefined;
}

/** The relative strand of an end pair: 1 straight, −1 reversed (the twist), 0 when either is unknown. */
function relative(a: Int8Array | undefined, b: Int8Array | undefined, i: number): number {
  const s = (a?.[i] ?? 0) * (b?.[i] ?? 0);
  return s > 0 ? 1 : s < 0 ? -1 : 0;
}

/**
 * Items to bins of `level` (`binPairs`): each bin pair with at least one pair, its count, by relative strand (when the
 * list has strands) and by `type` (`category`), a numeric field summed (`field`); bins in first-seen order.
 */
export function binPairs(t: PairsList, level: number, options: BinPairsOptions = {}): PairsBins {
  return binPairsRange(t, level, options, 0, t.aStart.length, new Binner()).finish();
}

/**
 * A binning that runs a chunk of items at a time (`step`), so a pyramid is built in time slices on the main thread;
 * `finish` gives the bins once every item is in.
 */
export function binPairsRange(
  t: PairsList,
  level: number,
  options: BinPairsOptions,
  from: number,
  to: number,
  binner: Binner,
) {
  const sysDict: string[] = [];
  const sysIndex = new Map<string, number>();
  const aSys = remap(t.aSystem, sysDict, sysIndex);
  const bSys = remap(t.bSystem, sysDict, sysIndex);
  const own = !aSys && options.system !== undefined && (options.systems || false);
  if (own) {
    sysDict.push(options.system!);
    sysIndex.set(options.system!, 0);
  }
  const catDict = options.category && t.type ? t.type.dict : undefined;
  const field = numericField(t, options.field);
  const strands = !!(t.aStrand || t.bStrand);
  const run = (until: number) => {
    for (let i = from; i < until; i++) {
      const as = aSys ? aSys[i] : 0;
      const bs = bSys ? bSys[i] : 0;
      const cat = catDict ? t.type!.codes[i] : -1;
      const ab = Math.floor(t.aStart[i] / level);
      const bb = Math.floor(t.bStart[i] / level);
      const s = field ? field[i] : 0;
      binner.add(as, bs, ab, bb, strands ? relative(t.aStrand, t.bStrand, i) : 0, cat, 1, Number.isNaN(s) ? 0 : s);
    }
    from = until;
  };
  return {
    /** Bins the items up to `until` (exclusive); `true` when every item is in. */
    step(until: number): boolean {
      run(Math.min(until, to));
      return from >= to;
    },
    finish(): PairsBins {
      run(to);
      return binner.finish(aSys || own ? sysDict : undefined, catDict, strands, field !== undefined);
    },
  };
}

/** A fresh accumulator for `binPairsRange`. */
export function pairsBinner(): Binner {
  return new Binner();
}

/**
 * Bins to bins `factor` times wider (`coarsenPairBins`): keyed as bins (both systems, both bins, strand, category),
 * counts and sums summed; exact, levels being nested.
 */
export function coarsenPairBins(t: PairsBins, factor: number): PairsBins {
  const b = new Binner();
  const aSys = t.aSystem?.codes;
  const bSys = t.bSystem?.codes;
  const cat = t.category?.codes;
  for (let i = 0; i < t.count.length; i++)
    b.add(
      aSys ? aSys[i] : 0,
      bSys ? bSys[i] : 0,
      Math.floor(t.aBin[i] / factor),
      Math.floor(t.bBin[i] / factor),
      t.strand ? t.strand[i] : 0,
      cat ? cat[i] : -1,
      t.count[i],
      t.sum ? t.sum[i] : 0,
    );
  // the systems' dictionaries are the input's: both ends' columns share one when the pyramid wrote them
  const out = b.finish(undefined, t.category?.dict, !!t.strand, !!t.sum);
  if (!t.aSystem || !t.bSystem) return out;
  return {
    ...out,
    aSystem: categoryFromCodes(t.aSystem.dict, b.aSys.subarray(0, b.n)),
    bSystem: categoryFromCodes(t.bSystem.dict, b.bSys.subarray(0, b.n)),
  };
}
