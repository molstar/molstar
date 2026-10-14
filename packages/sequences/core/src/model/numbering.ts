/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Numberings (seqspec:model#numberings): labels for the positions of a molecule or a unit; positions are indices.

import { CustomData, Issue, type NameId, type SeqIndex, type StringTable } from './common.js';
import { SortedArray } from './ordered-set.js';

/** An unlabeled position in a `labels` numbering. */
export const NO_LABEL = -2147483648;

export type NumberingMethod = 'offset' | 'segments' | 'labels' | 'names';
export interface NumberingBase {
  readonly name: NameId;
  readonly method: NumberingMethod;
  readonly isDefault: boolean;
  readonly custom: CustomData;
  /** Positions numbered: its molecule's length (an addition to seqspec:model#indices-and-labels: bounds `offset`). */
  readonly length: number;
}
/** label = index + first; with skipZero, labels jump from -1 to 1. */
export interface OffsetNumbering extends NumberingBase {
  readonly method: 'offset';
  readonly first: number;
  readonly skipZero: boolean;
}
/** Piecewise linear (UniProt by segments, SIFTS, isoforms): sorted by indexStart, may have holes. */
export interface SegmentNumbering extends NumberingBase {
  readonly method: 'segments';
  readonly indexStart: Int32Array;
  readonly labelStart: Int32Array;
  /** The spec's `length` per segment, renamed: `length` is NumberingBase's. */
  readonly segmentLength: Int32Array;
  /** Segments in labelStart order, built once. */
  readonly byLabel: Int32Array;
}
/** Integer labels with sparse suffixes (auth, Kabat, IMGT): "27A", "111.1". */
export interface LabelNumbering extends NumberingBase {
  readonly method: 'labels';
  /** NO_LABEL where unlabeled. */
  readonly num: Int32Array;
  /** Positions with a suffix, and its NameId ("A", ".1"). */
  readonly sufAt: SortedArray<SeqIndex>;
  readonly suffix: Int32Array;
  /** Each suffix's order (`rankOf`), read once at build: placing an absent endpoint makes no string (CM-09). */
  readonly sufRank: Float64Array;
  /** Built once: num → first index, dense from `lookupMin` (`lookupKeys` null) or sorted pairs. */
  readonly lookupMin: number;
  readonly lookupKeys: Int32Array | null;
  readonly lookupFirst: Int32Array;
  /** The next position with the same num, or -1. */
  readonly nextSame: Int32Array;
  /** Its lineage's table, to order suffixes in `rangeOf` (an addition to the spec). */
  readonly strings: StringTable;
}
/** A string per position: atom names, named pocket positions. */
export interface NameNumbering extends NumberingBase {
  readonly method: 'names';
  /** NameId per position, -1: unnamed. The spec's `name`, renamed: `name` is NumberingBase's. */
  readonly names: Int32Array;
  readonly byName: Map<number, number>;
}
export type Numbering = OffsetNumbering | SegmentNumbering | LabelNumbering | NameNumbering;

/** Reused by the caller. `num` NO_LABEL and `suffix` -1 for a name; `name` -1 for a number. */
export interface PositionLabel {
  num: number;
  /** The suffix's string id; `UNKNOWN_SUFFIX`: a written suffix no string of the set holds, ordered by `rank`. */
  suffix: NameId;
  name: NameId;
  /** With `UNKNOWN_SUFFIX`, the suffix's order (`Numbering.suffixOrder`): a reader never interns what it is asked. */
  rank?: number;
}

/**
 * A written suffix the set's strings do not hold (CM-11: a query reads its set, never interns into it): no position has
 * it, so the label is placed by order, the suffix's `rank`.
 */
export const UNKNOWN_SUFFIX = -2 as NameId;

export type NumberingInput = { name: string; isDefault?: boolean; custom?: object } & (
  | { method: 'offset'; first: number; skipZero?: boolean }
  | { method: 'segments'; segments: [number, number, number][] }
  | { method: 'labels'; labels: (number | string | null)[] }
  | { method: 'names'; names: (string | null)[] }
);

export interface NumberingSet {
  readonly count: number;
  readonly items: readonly Numbering[];
  byName(name: NameId): Numbering | null;
}

class NumberingList implements NumberingSet {
  readonly items: readonly Numbering[];
  constructor(items: readonly Numbering[]) {
    this.items = items;
  }
  get count() {
    return this.items.length;
  }
  byName(name: NameId): Numbering | null {
    for (let i = 0; i < this.items.length; i++) if (this.items[i].name === name) return this.items[i];
    return null;
  }
}

/**
 * A written label (seqspec:format#positions): an integer and an optional suffix, letters (an insertion code) or `.n`;
 * the model's one copy (CM-16: numberings, written annotation positions and queries read it). SeqQL's `sel` has its own,
 * as `query` imports nothing.
 */
export const LABEL = /^(-?[0-9]+)([A-Za-z]+|\.[0-9]+)?$/;
const RESERVED = new Set(['ordinal', 'row']);

function build(input: NumberingInput, length: number, strings: StringTable, issues: Issue[], path: string) {
  const err = (m: string) => {
    issues.push(Issue.error('model.numbering.invalid', path, `Numbering "${input.name}": ${m}`));
    return null;
  };
  const base = {
    name: strings.intern(input.name),
    isDefault: input.isDefault ?? false,
    custom: CustomData.create(input.custom),
    length,
  };
  if (input.method === 'offset') {
    if (!Number.isInteger(input.first)) return err('`first` must be an integer.');
    return { ...base, method: 'offset', first: input.first, skipZero: input.skipZero ?? false } as OffsetNumbering;
  }
  if (input.method === 'segments') {
    const n = input.segments.length;
    const indexStart = new Int32Array(n);
    const labelStart = new Int32Array(n);
    const segmentLength = new Int32Array(n);
    for (let k = 0; k < n; k++) {
      const [s, l, len] = input.segments[k];
      const prevStop = k ? indexStart[k - 1] + segmentLength[k - 1] : 0;
      if (len <= 0 || s < prevStop || s + len > length) {
        return err(`segment ${k} must have a positive length, follow the previous one and fit ${length} positions.`);
      }
      indexStart[k] = s;
      labelStart[k] = l;
      segmentLength[k] = len;
    }
    const byLabel = Int32Array.from(indexStart.keys()).sort((a, b) => labelStart[a] - labelStart[b]);
    for (let k = 1; k < n; k++) {
      const a = byLabel[k - 1];
      if (labelStart[byLabel[k]] < labelStart[a] + segmentLength[a]) {
        return err(`segment ${byLabel[k]} overlaps segment ${a} in label space; label ranges must be disjoint.`);
      }
    }
    return { ...base, method: 'segments', indexStart, labelStart, segmentLength, byLabel } as SegmentNumbering;
  }
  const values = input.method === 'labels' ? input.labels : input.names;
  if (values.length !== length) return err(`it has ${values.length} labels for ${length} positions.`);
  if (input.method === 'names') {
    const names = new Int32Array(length).fill(-1);
    const byName = new Map<number, number>();
    for (let i = 0; i < length; i++) {
      const v = input.names[i];
      if (v == null) continue;
      names[i] = strings.intern(v);
      if (!byName.has(names[i])) byName.set(names[i], i);
    }
    return { ...base, method: 'names', names, byName } as NameNumbering;
  }
  const num = new Int32Array(length).fill(NO_LABEL);
  const sufAt: number[] = [];
  const suffix: number[] = [];
  let lo = Infinity;
  let hi = -Infinity;
  for (let i = 0; i < length; i++) {
    const v = input.labels[i];
    if (v == null) continue;
    const m = typeof v === 'number' ? null : LABEL.exec(v);
    if (typeof v === 'number' ? !Number.isInteger(v) : !m) {
      return err(`label ${JSON.stringify(v)} at ${i} is not an integer with an optional suffix ("27A", "111.1").`);
    }
    const value = m ? Number(m[1]) : (v as number);
    // Held in an Int32Array, whose least value is NO_LABEL: a label outside is said, not wrapped (CM-14).
    if (!(value > NO_LABEL && value <= 0x7fffffff))
      return err(
        `label ${JSON.stringify(v)} at ${i} is outside the labels a numbering holds (-2147483647 to 2147483647).`,
      );
    num[i] = value;
    if (m?.[2]) {
      sufAt.push(i);
      suffix.push(strings.intern(m[2]));
    }
    lo = Math.min(lo, num[i]);
    hi = Math.max(hi, num[i]);
  }
  const dense = hi - lo + 1 <= 4 * length + 64;
  const keys = dense ? null : SortedArray.ofUnsortedArray(Int32Array.from(num).filter((x) => x !== NO_LABEL));
  const lookupFirst = new Int32Array(dense ? Math.max(0, hi - lo + 1) : (keys as Int32Array).length).fill(-1);
  const nextSame = new Int32Array(length).fill(-1);
  const last = new Int32Array(lookupFirst.length).fill(-1);
  for (let i = 0; i < length; i++) {
    if (num[i] === NO_LABEL) continue;
    const k = keys ? SortedArray.indexOf(keys, num[i]) : num[i] - lo;
    if (lookupFirst[k] < 0) lookupFirst[k] = i;
    else nextSame[last[k]] = i;
    last[k] = i;
  }
  return {
    ...base,
    method: 'labels',
    num,
    sufAt: Int32Array.from(sufAt),
    suffix: Int32Array.from(suffix),
    sufRank: Float64Array.from(suffix, (id) => rankOf(strings.get(id))),
    lookupMin: lo,
    lookupKeys: keys,
    lookupFirst,
    nextSame,
    strings,
  } as LabelNumbering;
}

/** A numbering set from inputs, `ordinal` built in when `withOrdinal`; a bad numbering is left out with an error. */
function buildSet(
  inputs: readonly NumberingInput[] | undefined,
  length: number,
  strings: StringTable,
  issues: Issue[],
  path: string,
  withOrdinal: boolean,
): NumberingSet {
  const items: Numbering[] = [];
  if (withOrdinal)
    items.push(build({ name: 'ordinal', method: 'offset', first: 1 }, length, strings, issues, path) as Numbering);
  for (let k = 0; k < (inputs?.length ?? 0); k++) {
    const input = (inputs as NumberingInput[])[k];
    const p = `${path}/${k}`;
    if (RESERVED.has(input.name) || items.some((n) => n.name === strings.find(input.name))) {
      issues.push(Issue.error('model.numbering.name', p, `Numbering name "${input.name}" is reserved or repeated.`));
      continue;
    }
    const n = build(input, length, strings, issues, p);
    if (n) items.push(n);
  }
  return new NumberingList(items);
}

function suffixAt(n: LabelNumbering, i: number): number {
  const k = SortedArray.indexOf(n.sufAt, i);
  return k < 0 ? -1 : n.suffix[k];
}

function indexOf(n: Numbering, num: number, suffix: NameId): SeqIndex {
  switch (n.method) {
    case 'offset': {
      if (suffix !== -1 || (n.skipZero && n.first <= 0 && num === 0)) return -1 as SeqIndex;
      const i = num - n.first - (n.skipZero && n.first <= 0 && num > 0 ? 1 : 0);
      return (i >= 0 && i < n.length ? i : -1) as SeqIndex;
    }
    case 'segments': {
      if (suffix !== -1) return -1 as SeqIndex;
      // The last segment (in label order) starting at or before num.
      let lo = 0;
      let hi = n.byLabel.length;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (n.labelStart[n.byLabel[mid]] <= num) lo = mid + 1;
        else hi = mid;
      }
      if (lo === 0) return -1 as SeqIndex;
      const s = n.byLabel[lo - 1];
      const d = num - n.labelStart[s];
      return (d < n.segmentLength[s] ? n.indexStart[s] + d : -1) as SeqIndex;
    }
    case 'labels': {
      const k = n.lookupKeys ? SortedArray.indexOf(n.lookupKeys, num) : num - n.lookupMin;
      let i = k >= 0 && k < n.lookupFirst.length ? n.lookupFirst[k] : -1;
      while (i >= 0 && suffixAt(n, i) !== suffix) i = n.nextSame[i];
      return i as SeqIndex;
    }
    default:
      return -1 as SeqIndex;
  }
}

function indexOfName(n: NameNumbering, name: NameId): SeqIndex {
  return (n.byName.get(name) ?? -1) as SeqIndex;
}

function labelOf(n: Numbering, i: SeqIndex, out: PositionLabel): boolean {
  out.num = NO_LABEL;
  out.suffix = -1 as NameId;
  out.name = -1 as NameId;
  if (i < 0 || i >= n.length) return false;
  switch (n.method) {
    case 'offset': {
      const v = i + n.first;
      out.num = n.skipZero && n.first <= 0 && v >= 0 ? v + 1 : v;
      return true;
    }
    case 'segments': {
      const k = SortedArray.lowerBound(n.indexStart, i + 1) - 1;
      if (k < 0 || i >= n.indexStart[k] + n.segmentLength[k]) return false;
      out.num = n.labelStart[k] + i - n.indexStart[k];
      return true;
    }
    case 'labels':
      if (n.num[i] === NO_LABEL) return false;
      out.num = n.num[i];
      out.suffix = suffixAt(n, i) as NameId;
      return true;
    case 'names':
      out.name = n.names[i] as NameId;
      return out.name >= 0;
  }
}

/**
 * Suffix order for placing absent endpoints: none, then letters (by length, then alphabetically: A … Z, AA … ZZ, as
 * insertion codes run out past Z), then `.n` by number.
 */
function suffixRank(n: Numbering, l: PositionLabel): number {
  if (l.suffix === UNKNOWN_SUFFIX) return l.rank ?? 1;
  if (l.suffix < 0) return 0;
  if (n.method !== 'labels') return 1;
  return rankOf(n.strings.get(l.suffix));
}
/** A suffix's place in that order, from its text. */
function rankOf(s: string): number {
  if (s[0] === '.') return 1e12 + Number(s.slice(1));
  const u = s.toUpperCase();
  let r = 0;
  for (let k = 0; k < u.length && k < 6; k++) r = r * 27 + (u.charCodeAt(k) - 64);
  return r;
}

const scratch: PositionLabel = { num: 0, suffix: -1 as NameId, name: -1 as NameId };

/**
 * A label range's ends placed by seqspec:format#positions, their order not judged: `out[0]` the index of `beg` (or,
 * absent, the first position at or after it), `out[1]` of `end` (or the last at or before it); -1 where none is
 * placed. Tells a range written backwards (both placed, the end first) from a label that is not there (CF-01).
 */
function placeRange(n: Numbering, beg: PositionLabel, end: PositionLabel, out: Int32Array): void {
  let start: number;
  let last: number;
  if (n.method === 'names') {
    start = indexOfName(n, beg.name);
    last = indexOfName(n, end.name);
  } else {
    start = indexOf(n, beg.num, beg.suffix);
    last = indexOf(n, end.num, end.suffix);
    if (start < 0 || last < 0) {
      const begRank = suffixRank(n, beg);
      const endRank = suffixRank(n, end);
      const findStart = start < 0;
      const findLast = last < 0;
      // The positions' suffix ranks from the build's table, in step with them: no string per position (CM-09).
      const sufAt = n.method === 'labels' ? n.sufAt : null;
      let j = 0;
      for (let i = 0; i < n.length; i++) {
        if (!labelOf(n, i as SeqIndex, scratch)) continue;
        while (sufAt && j < sufAt.length && sufAt[j] < i) j++;
        const rank = sufAt && j < sufAt.length && sufAt[j] === i ? (n as LabelNumbering).sufRank[j] : 0;
        if (findStart && start < 0 && (scratch.num - beg.num || rank - begRank) >= 0) start = i;
        if (findLast && (scratch.num - end.num || rank - endRank) <= 0) last = i;
      }
    }
  }
  out[0] = start;
  out[1] = last;
}

const ends = new Int32Array(2);

/**
 * Whether `beg` is written after `end` (CF-01: a range written backwards): by index when both are labels of the
 * numbering (seqspec:model#numberings: IMGT's 112.1 comes before its 112, a circular permutant's 53 before its 1), by
 * (num, suffix) as `placeRange` orders an absent one. False when either is not a label of the numbering's kind.
 */
function backwards(n: Numbering, beg: PositionLabel, end: PositionLabel): boolean {
  if (n.method === 'names') {
    const [a, b] = [indexOfName(n, beg.name), indexOfName(n, end.name)];
    return a >= 0 && b >= 0 && b < a;
  }
  if (beg.num === NO_LABEL || end.num === NO_LABEL) return false;
  const [a, b] = [indexOf(n, beg.num, beg.suffix), indexOf(n, end.num, end.suffix)];
  if (a >= 0 && b >= 0) return b < a;
  return (beg.num - end.num || suffixRank(n, beg) - suffixRank(n, end)) > 0;
}

/**
 * A label range by seqspec:format#positions: index order between present endpoints; an absent endpoint is placed
 * by (num, suffix): `beg` at the first position at or after it, `end` at the last at or before it.
 * Writes [start, stop) to `out`; false when it covers nothing. `cyclic`: `beg` after `end` wraps the origin (start > stop).
 */
function rangeOf(n: Numbering, beg: PositionLabel, end: PositionLabel, out: Int32Array, cyclic = false): boolean {
  placeRange(n, beg, end, ends);
  const start = ends[0];
  const last = ends[1];
  if (start < 0 || last < 0 || (last < start && !(cyclic && last + 1 < start))) return false;
  out[0] = start;
  out[1] = last + 1;
  return true;
}

function format(n: Numbering, i: SeqIndex, strings: StringTable): string {
  if (!labelOf(n, i, scratch)) return '';
  return n.method === 'names' ? strings.get(scratch.name) : `${scratch.num}${strings.get(scratch.suffix)}`;
}

/** A set of built numberings (row spaces' `ordinal` and `row`). */
const setOf = (items: readonly Numbering[]): NumberingSet => new NumberingList(items);

export const Numbering = {
  build,
  buildSet,
  setOf,
  indexOf,
  indexOfName,
  labelOf,
  placeRange,
  backwards,
  rangeOf,
  format,
  /** A suffix's order among labels' suffixes (`UNKNOWN_SUFFIX`'s `rank`). */
  suffixOrder: rankOf,
};
