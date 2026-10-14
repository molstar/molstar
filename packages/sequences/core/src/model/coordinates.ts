/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Coordinate systems and coordinate maps (seqspec:model#coordinate-maps; seqspec:model#coordinate-systems). A
// coordinate map is one segment table whatever made it (identity, a unit step, a row onto its columns, a mapping's
// method, a codon's CDS, a composed path): each segment pairs the from sub-positions [at, at + length), q =
// span·element + phase, with `length` to positions from `toLo`, ascending (`dir` 1) or descending (-1, strand -). Built
// at the edges; lookups and range walks read the arrays and allocate nothing (the callbacks take a context, as Mol*'s
// iterators do).

import type { Alignment } from "./alignment.ts";
import type { RowIndex, UnitId } from "./common.ts";
import type { SeqSet } from "./set.ts";

export const CoordinateSystemKind = {
  molecule: 0,
  unit: 1,
  row: 2,
  columns: 3,
  complex: 4,
  external: 5,
  item: 6,
} as const;
/** id: a MoleculeId, UnitId, AlignmentId, ComplexId, AnnotationId or registry id; row: a row's or an item's, else -1. */
export interface CoordinateSystem {
  kind: number;
  id: number;
  row: number;
}
export interface Coordinate {
  system: CoordinateSystem;
  index: number;
}

export interface CoordinateMap {
  readonly from: CoordinateSystem;
  readonly to: CoordinateSystem;
  /** To positions per from element: 1, or 3 through a codon step. */
  readonly span: number;
  /** Segments, sorted by `at`. */
  readonly count: number;
  /** The segment's first from sub-position, q = span·element + phase. */
  readonly at: Int32Array;
  /** Sub-positions in it, = to positions. */
  readonly length: Int32Array;
  /** The lowest to position it covers. */
  readonly toLo: Int32Array;
  /** 1: to ascends with q; -1: descends (strand -). */
  readonly dir: Int8Array;
  /** Segment order by `toLo` (then `at`): segments over one to range sit together, copies in from order. */
  readonly byTo: Int32Array;
  /** Per `byTo` entry: the highest to stop of it and those before, so a to → from search stops where none reaches. */
  readonly byToMax: Int32Array;
  /**
   * 1: `at` order is to order (to ranges ascending and disjoint), -1: descending; 0: neither, a range walk reads `byTo`;
   * 2: from ranges overlap too (a composed path with a source used twice at both ends; `toAt` gives the first).
   */
  readonly order: number;
  /** Support hulls, in elements and in to positions. */
  readonly fromStart: number;
  readonly fromStop: number;
  readonly toStart: number;
  readonly toStop: number;
}
/** Reused by the caller. */
export interface MapHit {
  element: number;
  phase: number;
}

type Run<C> = (start: number, stop: number, ctx: C) => void;

/** A map from segments as (at, length, toLo, dir) quadruples in any order: sorted by `at`, collinear neighbors merged. */
function create(from: CoordinateSystem, to: CoordinateSystem, span: number, quads: readonly number[]): CoordinateMap {
  const order = Array.from({ length: quads.length >> 2 }, (_, i) => 4 * i)
    .filter((i) => quads[i + 1] > 0)
    .sort((a, b) => quads[a] - quads[b] || quads[a + 2] - quads[b + 2]);
  const s: number[][] = [[], [], [], []];
  for (const i of order) {
    const [a, l, t, d] = [quads[i], quads[i + 1], quads[i + 2], quads[i + 3]];
    const k = s[0].length - 1;
    if (k >= 0 && s[3][k] === d && s[0][k] + s[1][k] === a && (d > 0 ? s[2][k] + s[1][k] === t : t + l === s[2][k])) {
      s[1][k] += l;
      if (d < 0) s[2][k] = t;
    } else for (const [j, x] of [a, l, t, d].entries()) s[j].push(x);
  }
  const [at, length, toLo, dir] = [
    Int32Array.from(s[0]),
    Int32Array.from(s[1]),
    Int32Array.from(s[2]),
    Int8Array.from(s[3]),
  ];
  const count = at.length;
  const byTo = Int32Array.from(at.keys()).sort((a, b) => toLo[a] - toLo[b] || at[a] - at[b]);
  const byToMax = new Int32Array(count);
  let [up, down, overlap, fromStop, toStart, toStop] = [true, true, false, 0, count ? 0x7fffffff : 0, 0];
  for (let k = 0; k < count; k++) {
    byToMax[k] = Math.max(k ? byToMax[k - 1] : 0, toLo[byTo[k]] + length[byTo[k]]);
    if (k && at[k] < fromStop) overlap = true;
    if (k && toLo[k] < toLo[k - 1] + length[k - 1]) up = false;
    if (k && toLo[k] + length[k] > toLo[k - 1]) down = false;
    fromStop = Math.max(fromStop, at[k] + length[k]);
    toStart = Math.min(toStart, toLo[k]);
    toStop = Math.max(toStop, toLo[k] + length[k]);
  }
  const fromStart = count ? Math.floor(at[0] / span) : 0;
  fromStop = Math.ceil(fromStop / span);
  const o = overlap ? 2 : up ? 1 : down ? -1 : 0;
  return {
    from,
    to,
    span,
    count,
    at,
    length,
    toLo,
    dir,
    byTo,
    byToMax,
    order: o,
    fromStart,
    fromStop,
    toStart,
    toStop,
  };
}

/** The identity over a domain, [start, stop). */
function identity(system: CoordinateSystem, [start, stop]: readonly [number, number]): CoordinateMap {
  return create(system, system, 1, [start, stop - start, start, 1]);
}

/** The same relation read to → from (span 1 only: a span map keeps its protein side as `from`). Edges. */
function invert(m: CoordinateMap): CoordinateMap {
  if (m.span !== 1) throw new Error("Coords.invert: a span map keeps its protein side as `from`.");
  const quads: number[] = [];
  for (let k = 0; k < m.count; k++) quads.push(m.toLo[k], m.length[k], m.at[k], m.dir[k]);
  return create(m.to, m.from, 1, quads);
}

/**
 * `a` then `b` (a's to is b's from), as a map from a's `from` to `to`, clipped to `domain` there. At most one of them
 * has span 3; a span-1 map before a codon step is read per sub-position. Sorted sweep: O((m + n) log n). Edges.
 */
function compose(
  a: CoordinateMap,
  b: CoordinateMap,
  to: CoordinateSystem,
  [lo, hi]: readonly [number, number],
): CoordinateMap {
  if (a.span > 1 && b.span > 1) throw new Error("Coords.compose: a path holds at most one codon step.");
  const x = b.span > a.span ? b.span : 1;
  const quads: number[] = [];
  for (let i = 0; i < a.count; i++) {
    if (x > 1 && a.dir[i] < 0 && a.length[i] > 1)
      throw new Error("Coords.compose: a strand - step before a codon step.");
    // Read per sub-position, one element keeps its phases in order whatever its direction.
    const [aAt, y0, y1, ad] = [x * a.at[i], x * a.toLo[i], x * (a.toLo[i] + a.length[i]), x > 1 ? 1 : a.dir[i]];
    for (let j = lastBelow(b.at, b.count, y1); j >= 0; j--) {
      const [bs, bt] = [b.at[j], b.at[j] + b.length[j]];
      if (bt <= y0) {
        if (b.order === 2) continue;
        break;
      }
      const [u0, u1] = [Math.max(y0, bs), Math.min(y1, bt)];
      if (u0 >= u1) continue;
      const z0 = b.dir[j] > 0 ? b.toLo[j] + (u0 - bs) : b.toLo[j] + b.length[j] - (u1 - bs);
      const [c0, c1] = [Math.max(z0, lo), Math.min(z0 + u1 - u0, hi)];
      if (c0 >= c1) continue;
      const d = ad * b.dir[j];
      const f0 = ad > 0 ? aAt + (u0 - y0) : aAt + (y1 - u1);
      quads.push(d > 0 ? f0 + (c0 - z0) : f0 + (z0 + u1 - u0 - c1), c1 - c0, c0, d);
    }
  }
  return create(a.from, to, a.span * b.span, quads);
}

/** The last index k < count with sorted[k] < x; -1 when none. */
function lastBelow(sorted: Int32Array, count: number, x: number): number {
  let lo = 0;
  let hi = count;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (sorted[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo - 1;
}

/** The last `byTo` entry whose segment starts at or below `to`; -1 when none. */
function lastByTo(m: CoordinateMap, to: number): number {
  let lo = 0;
  let hi = m.count;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (m.toLo[m.byTo[mid]] <= to) lo = mid + 1;
    else hi = mid;
  }
  return lo - 1;
}

/** The segment holding sub-position q, -1: none. */
function segmentAt(m: CoordinateMap, q: number): number {
  for (let k = lastBelow(m.at, m.count, q + 1); k >= 0; k--) {
    if (q < m.at[k] + m.length[k]) return k;
    if (m.order !== 2) break;
  }
  return -1;
}

/** The to position of sub-position q in segment k. */
const toOf = (m: CoordinateMap, k: number, q: number) =>
  m.dir[k] > 0 ? m.toLo[k] + (q - m.at[k]) : m.toLo[k] + m.length[k] - 1 - (q - m.at[k]);

/** The to position of an element's phase; -1: none. */
function toAt(m: CoordinateMap, element: number, phase: number): number {
  const q = m.span * element + phase;
  const k = segmentAt(m, q);
  return k < 0 ? -1 : toOf(m, k, q);
}

/** Every from element holding a to position (two for a source used twice), in from order: the count visited. */
function forEachFromAt<C>(m: CoordinateMap, to: number, out: MapHit, visit: (h: MapHit, c: C) => void, c: C): number {
  const hi = lastByTo(m, to);
  if (hi < 0) return 0;
  let lo = hi;
  while (lo > 0 && m.byToMax[lo - 1] > to) lo--;
  let n = 0;
  for (let k = lo; k <= hi; k++) {
    const s = m.byTo[k];
    if (to >= m.toLo[s] + m.length[s]) continue;
    const q = m.at[s] + (m.dir[s] > 0 ? to - m.toLo[s] : m.toLo[s] + m.length[s] - 1 - to);
    out.element = (q / m.span) | 0;
    out.phase = q - out.element * m.span;
    visit(out, c);
    n++;
  }
  return n;
}

/** The pending run of each range walk in progress, by depth (a visitor may walk another map): touching runs merge. */
const runs = new Int32Array(64).fill(-1);
let depth = 0;
function put<C>(d: number, s: number, t: number, visit: Run<C>, c: C) {
  const i = 2 * d;
  const j = i + 1;
  if (runs[i] >= 0 && s <= runs[j] && t >= runs[i]) {
    runs[i] = Math.min(runs[i], s);
    runs[j] = Math.max(runs[j], t);
    return;
  }
  if (runs[i] >= 0) visit(runs[i], runs[j], c);
  runs[i] = s;
  runs[j] = t;
}
function flush<C>(d: number, visit: Run<C>, c: C) {
  const s = runs[2 * d];
  runs[2 * d] = -1;
  if (s >= 0) visit(s, runs[2 * d + 1], c);
}
/** Segment k's sub-positions [q0, q1) clipped to it, as a to run. */
function putTo<C>(d: number, m: CoordinateMap, k: number, q0: number, q1: number, visit: Run<C>, c: C) {
  const a = Math.max(q0, m.at[k]);
  const b = Math.min(q1, m.at[k] + m.length[k]);
  if (a >= b) return;
  if (m.dir[k] > 0) put(d, toOf(m, k, a), toOf(m, k, b - 1) + 1, visit, c);
  else put(d, toOf(m, k, b - 1), toOf(m, k, a) + 1, visit, c);
}

/** Per nesting depth, the clipped runs of an order-0 or -2 map, sorted before they are put (grown, never shrunk). */
const clipped: Int32Array[] = [];

/** Segment s clipped to [q0, q1) as a to run into `out` at pair `n`; false when it holds none. */
function clipTo(m: CoordinateMap, s: number, q0: number, q1: number, out: Int32Array, n: number): boolean {
  const a = Math.max(q0, m.at[s]);
  const b = Math.min(q1, m.at[s] + m.length[s]);
  if (a >= b) return false;
  out[2 * n] = m.dir[s] > 0 ? toOf(m, s, a) : toOf(m, s, b - 1);
  out[2 * n + 1] = m.dir[s] > 0 ? toOf(m, s, b - 1) + 1 : toOf(m, s, a) + 1;
  return true;
}

/**
 * From elements [start, stop) as to runs of consecutive positions, in to order (a codon split by an intron: two runs).
 * A boundary (start === stop) before element k is visited as (b, b): before k's first to position in reading order (on
 * strand -, the one after it by index); k unmapped or at the end, after the last position of k − 1.
 */
function forEachTo<C>(m: CoordinateMap, start: number, stop: number, visit: Run<C>, c: C): void {
  if (start === stop) {
    boundary(m, start, visit, c);
    return;
  }
  const q0 = m.span * start;
  const q1 = m.span * stop;
  const d = depth++;
  try {
    if (m.order === 0 || m.order === 2) {
      // Segments in to order clip to runs out of it (a tandem map: [90, 100) then [0, 10)): sorted by start, then
      // put, which merges what overlaps or touches (CM-08). Insertion sort in a per-depth scratch: nothing allocated
      // once it has grown.
      let out = clipped[d];
      if (!out || out.length < 2 * m.count) {
        out = new Int32Array(Math.max(16, 2 * m.count));
        clipped[d] = out;
      }
      let n = 0;
      for (let k = 0; k < m.count; k++) if (clipTo(m, m.byTo[k], q0, q1, out, n)) n++;
      for (let i = 1; i < n; i++) {
        const s = out[2 * i];
        const t = out[2 * i + 1];
        let j = i - 1;
        for (; j >= 0 && out[2 * j] > s; j--) {
          out[2 * j + 2] = out[2 * j];
          out[2 * j + 3] = out[2 * j + 1];
        }
        out[2 * j + 2] = s;
        out[2 * j + 3] = t;
      }
      for (let i = 0; i < n; i++) put(d, out[2 * i], out[2 * i + 1], visit, c);
    } else {
      const last = lastBelow(m.at, m.count, q1);
      let first = Math.max(0, last);
      while (first > 0 && m.at[first - 1] + m.length[first - 1] > q0) first--;
      if (m.order > 0) for (let k = first; k <= last; k++) putTo(d, m, k, q0, q1, visit, c);
      else for (let k = last; k >= first; k--) putTo(d, m, k, q0, q1, visit, c);
    }
    flush(d, visit, c);
  } finally {
    runs[2 * d] = -1;
    depth = d;
  }
}

function boundary<C>(m: CoordinateMap, k: number, visit: Run<C>, c: C) {
  for (let p = 0; p < m.span; p++) {
    const s = segmentAt(m, m.span * k + p);
    if (s < 0) continue;
    const t = toOf(m, s, m.span * k + p);
    visit(m.dir[s] > 0 ? t : t + 1, m.dir[s] > 0 ? t : t + 1, c);
    return;
  }
  for (let p = m.span - 1; p >= 0 && k > 0; p--) {
    const s = segmentAt(m, m.span * (k - 1) + p);
    if (s < 0) continue;
    const t = toOf(m, s, m.span * (k - 1) + p);
    visit(m.dir[s] > 0 ? t + 1 : t, m.dir[s] > 0 ? t + 1 : t, c);
    return;
  }
}

/** To positions [start, stop) as runs of every from element holding them (the inverse image, all copies). */
function forEachFrom<C>(m: CoordinateMap, start: number, stop: number, visit: Run<C>, c: C): void {
  const hi = lastByTo(m, stop - 1);
  let lo = Math.max(0, hi);
  while (lo > 0 && m.byToMax[lo - 1] > start) lo--;
  const d = depth++;
  try {
    for (let k = lo; k <= hi && start < stop; k++) {
      const s = m.byTo[k];
      const x0 = Math.max(start, m.toLo[s]);
      const x1 = Math.min(stop, m.toLo[s] + m.length[s]);
      if (x0 >= x1) continue;
      const qa = m.at[s] + (m.dir[s] > 0 ? x0 - m.toLo[s] : m.toLo[s] + m.length[s] - x1);
      put(d, (qa / m.span) | 0, (((qa + x1 - x0 - 1) / m.span) | 0) + 1, visit, c);
    }
    flush(d, visit, c);
  } finally {
    runs[2 * d] = -1;
    depth = d;
  }
}

/**
 * An insert's next pair (seqspec:model#projection): for `[start, stop)` on one side of `m` (`from`: its elements, else
 * its `to` positions) between positions the map pairs, one before `start`, the first paired at `stop` or after; -1: not
 * between pairs. Binary searches over the segments (`at` order, `byTo` order); allocation-free.
 */
function pairedAfter(m: CoordinateMap, start: number, stop: number, from: boolean): number {
  if (m.count === 0) return -1;
  if (from) {
    if (Math.floor(m.at[0] / m.span) >= start) return -1;
    if (segmentAt(m, m.span * stop) >= 0) return stop;
    const k = lastBelow(m.at, m.count, m.span * stop) + 1;
    return k < m.count ? Math.floor(m.at[k] / m.span) : -1;
  }
  if (m.toLo[m.byTo[0]] >= start) return -1;
  const hi = lastByTo(m, stop);
  if (hi >= 0 && m.byToMax[hi] > stop) return stop;
  return hi + 1 < m.count ? m.toLo[m.byTo[hi + 1]] : -1;
}

/**
 * A boundary on `m`'s `to` side, before position `at`, as the boundary on its `from` side, visited as (k, k): before the
 * element holding `at` (after it on a reversed segment); `at` not held, after the element holding `at − 1` (before it
 * on a reversed segment); within a codon, after its residue; neither held, none (false). A position held twice: its
 * first segment's.
 */
function boundaryFrom<C>(m: CoordinateMap, at: number, visit: Run<C>, c: C): boolean {
  for (let p = 0; p < 2; p++) {
    const to = at - p;
    const hi = to < 0 ? -1 : lastByTo(m, to);
    let lo = hi;
    while (lo > 0 && m.byToMax[lo - 1] > to) lo--;
    for (let k = lo; k <= hi && hi >= 0; k++) {
      const s = m.byTo[k];
      if (to >= m.toLo[s] + m.length[s]) continue;
      const q = m.at[s] + (m.dir[s] > 0 ? to - m.toLo[s] : m.toLo[s] + m.length[s] - 1 - to);
      // The boundary in reading order, a sub-position; inside a codon, after its residue.
      const r = q + ((p === 0) === m.dir[s] > 0 ? 0 : 1);
      const b = Math.ceil(r / m.span);
      visit(b, b, c);
      return true;
    }
  }
  return false;
}

/** A unit's step to its molecule: the identity within its range. */
function ofUnit(set: SeqSet, u: UnitId): CoordinateMap {
  const [s, t] = [set.units.rangeStart[u], set.units.rangeStop[u]];
  const unit = { kind: CoordinateSystemKind.unit, id: u, row: -1 };
  return create(unit, { kind: CoordinateSystemKind.molecule, id: set.units.molecule[u], row: -1 }, 1, [s, t - s, s, 1]);
}

/** A row onto its alignment's columns: one segment per ungapped run, A3M inserts mapped to nothing; null: not held. */
function ofRow(a: Alignment, row: RowIndex): CoordinateMap | null {
  const cols = a.rows.posMap(row);
  if (!cols) return null;
  const quads: number[] = [];
  for (let p = 0; p < cols.length; p++) {
    if (cols[p] < 0) continue;
    const k = quads.length - 4;
    if (k >= 0 && quads[k] + quads[k + 1] === p && quads[k + 2] + quads[k + 1] === cols[p]) quads[k + 1]++;
    else quads.push(p, 1, cols[p], 1);
  }
  const columns = { kind: CoordinateSystemKind.columns, id: a.id, row: -1 };
  return create({ kind: CoordinateSystemKind.row, id: a.id, row }, columns, 1, quads);
}

// Cyclic coordinate systems (seqspec:model#cyclic-systems): the domain [0, L), positions modulo L. An extent never
// wraps: a range across the origin is written start > stop and held as its two runs, [start, L) then [0, stop).

/** The distance between two positions of a cyclic system of `length`: min(|a − b|, L − |a − b|). */
function cyclicDistance(a: number, b: number, length: number): number {
  const d = Math.abs(a - b) % length;
  return Math.min(d, length - d);
}

/**
 * A range of a cyclic system of `length` as runs into `out` (pairs), their count: [start, stop) as written, or across
 * the origin (start > stop) [start, L) then [0, stop), an empty run left out. Allocation-free.
 */
function cyclicRuns(start: number, stop: number, length: number, out: Int32Array): number {
  if (start <= stop) {
    out[0] = start;
    out[1] = stop;
    return 1;
  }
  out[0] = start;
  out[1] = length;
  if (stop <= 0) return 1;
  out[2] = 0;
  out[3] = stop;
  return 2;
}

/**
 * The shortest arc of a cyclic system of `length` holding the runs [s0, t0, s1, t1, …] (sorted by start, within
 * [0, L)): the circle less its largest gap, as one run, or two across the origin, into `out`; their count (0: none).
 * Ties keep the arc that does not wrap. Focus, scroll-to and a projection's extent read it, so a feature across the
 * origin is one short arc, never the whole system. Allocation-free, one pass.
 */
function cyclicHull(runs: Int32Array, count: number, length: number, out: Int32Array): number {
  if (count <= 0) return 0;
  // The gap across the origin, from the furthest stop round to the first start, is the one a linear hull leaves out.
  let reach = runs[1];
  for (let k = 1; k < count; k++) reach = Math.max(reach, runs[2 * k + 1]);
  let gap = runs[0] + length - reach;
  let s = runs[0];
  let t = reach;
  let end = runs[1];
  for (let k = 1; k < count; k++) {
    const g = runs[2 * k] - end;
    if (g > gap) {
      gap = g;
      s = runs[2 * k];
      t = end;
    }
    end = Math.max(end, runs[2 * k + 1]);
  }
  return cyclicRuns(s, t, length, out);
}

const sameSystem = (a: CoordinateSystem, b: CoordinateSystem) => a.kind === b.kind && a.id === b.id && a.row === b.row;

/** Two maps relating the same elements the same way (their canonical segment tables are equal). */
function areEqual(a: CoordinateMap, b: CoordinateMap): boolean {
  if (!sameSystem(a.from, b.from) || !sameSystem(a.to, b.to) || a.span !== b.span || a.count !== b.count) return false;
  for (let k = 0; k < a.count; k++) {
    if (a.at[k] !== b.at[k] || a.length[k] !== b.length[k] || a.toLo[k] !== b.toLo[k] || a.dir[k] !== b.dir[k]) {
      return false;
    }
  }
  return true;
}

export const Coords = {
  create,
  identity,
  invert,
  compose,
  toAt,
  forEachFromAt,
  forEachTo,
  forEachFrom,
  pairedAfter,
  boundaryFrom,
  ofUnit,
  ofRow,
  areEqual,
  sameSystem,
  cyclicDistance,
  cyclicRuns,
  cyclicHull,
};
