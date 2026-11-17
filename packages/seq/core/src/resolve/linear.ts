/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The one-to-one part of a map and its encodings (_spec/README.md "Resolve", Maps): implicit (arithmetic), runs
// (segments) and expanded (a position each, for a region). Each answers the same interface: the image of a range
// forward and back as runs, the whole map as runs when it is bounded, and its inverse. Composition composes encodings
// (two implicit maps give an implicit one, anything else gives runs) and never expands one.

/**
 * Source positions `[from, from + length)` mapped to the target: `from + i` is `to + dir * i`. Direction is part of
 * every run, so a reverse strand and an inversion are the same mechanism.
 */
export interface Run {
  readonly from: number;
  readonly to: number;
  readonly length: number;
  readonly dir: 1 | -1;
}

export type Encoding = 'implicit' | 'runs' | 'expanded';

/** A partial map between the positions of two systems, at most one image per source position. */
export interface Linear {
  readonly encoding: Encoding;
  /** The source positions it may map, `[lo, hi)`; `hi` is Infinity on an open system. */
  readonly lo: number;
  readonly hi: number;
  /** The images of `[start, end)` as runs, sorted by source; what no run holds maps nowhere. */
  forward(start: number, end: number): Run[];
  /** The preimages of target positions `[start, end)` as runs from the target to the source, sorted by target. */
  backward(start: number, end: number): Run[];
  /** The whole map as canonical runs, or undefined when it is unbounded (an offset on an open system). */
  runs(): readonly Run[] | undefined;
  inverse(): Linear;
}

/** The target positions of a run, `[lo, hi)`. */
export function runTargets(r: Run): [number, number] {
  return r.dir === 1 ? [r.to, r.to + r.length] : [r.to - r.length + 1, r.to + 1];
}

/** The run from the target back to the source. */
export function invertRun(r: Run): Run {
  if (r.dir === 1) return { from: r.to, to: r.from, length: r.length, dir: 1 };
  return { from: r.to - r.length + 1, to: r.from + r.length - 1, length: r.length, dir: -1 };
}

/** The part of a run whose source is in `[start, end)`, keeping its direction, or undefined. */
export function clipRun(r: Run, start: number, end: number): Run | undefined {
  const s = Math.max(r.from, start);
  const e = Math.min(r.from + r.length, end);
  if (e <= s) return undefined;
  return { from: s, to: r.to + r.dir * (s - r.from), length: e - s, dir: r.dir };
}

const byFrom = (a: Run, b: Run) => a.from - b.from || a.to - b.to;

/**
 * Runs in their one canonical form: sorted by source, and joined greedily position by position, so equal maps give
 * equal lists whatever runs they were built from. A run of one position takes its direction from a neighbor it joins;
 * runs over the same source positions (a map read backwards that is one-to-many) are kept apart.
 */
export function canonicalRuns(input: readonly Run[]): Run[] {
  const runs = input.filter((r) => r.length > 0).sort(byFrom);
  const out: Run[] = [];
  let cur: { from: number; to: number; length: number; dir: 1 | -1 } | undefined;
  const push = () => {
    if (cur) out.push({ from: cur.from, to: cur.to, length: cur.length, dir: cur.dir });
  };
  for (const r of runs) {
    if (cur && r.from === cur.from + cur.length) {
      const last = cur.to + cur.dir * (cur.length - 1);
      const step = r.to - last;
      if (step === 1 || step === -1) {
        const dir = cur.length === 1 ? (step as 1 | -1) : cur.dir;
        if (step === dir) {
          cur.dir = dir;
          if (r.length === 1 || r.dir === dir) {
            cur.length += r.length;
            continue;
          }
          // the run's first position joins, the rest turns
          cur.length += 1;
          push();
          cur = { from: r.from + 1, to: r.to + r.dir, length: r.length - 1, dir: r.dir };
          continue;
        }
      }
    }
    push();
    cur = { from: r.from, to: r.to, length: r.length, dir: r.dir };
  }
  push();
  return out;
}

/** Whether two canonical run lists are the same map (a run of one position has no direction). */
export function runsEqual(a: readonly Run[], b: readonly Run[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    const x = a[i];
    const y = b[i];
    if (x.from !== y.from || x.to !== y.to || x.length !== y.length) return false;
    if (x.length > 1 && x.dir !== y.dir) return false;
  }
  return true;
}

/** A run cut where it crosses the origin of a cyclic target of `size` positions (EC-11). */
export function wrapRun(r: Run, size: number): Run[] {
  const out: Run[] = [];
  let from = r.from;
  let left = r.length;
  let to = (((r.to % size) + size) % size) | 0;
  while (left > 0) {
    const room = r.dir === 1 ? size - to : to + 1;
    const n = Math.min(left, room);
    out.push({ from, to, length: n, dir: r.dir });
    from += n;
    left -= n;
    to = r.dir === 1 ? 0 : size - 1;
  }
  return out;
}

// ---- implicit ----

/** `k ↦ d + dir·k` on `[lo, hi)`, taken modulo `size` on a cyclic target. Nothing stored. */
export class ImplicitLinear implements Linear {
  readonly encoding = 'implicit';
  constructor(
    readonly d: number,
    readonly dir: 1 | -1,
    readonly lo: number,
    readonly hi: number,
    /** The target's length when it is cyclic: images wrap. Needs a bounded domain. */
    readonly size?: number,
  ) {
    if (size !== undefined && !Number.isFinite(hi - lo)) throw new Error('A wrapping map needs a bounded domain');
  }

  forward(start: number, end: number): Run[] {
    const s = Math.max(start, this.lo);
    const e = Math.min(end, this.hi);
    if (e <= s) return [];
    const run: Run = { from: s, to: this.d + this.dir * s, length: e - s, dir: this.dir };
    return this.size === undefined ? [run] : wrapRun(run, this.size);
  }

  backward(start: number, end: number): Run[] {
    return this.inverse().forward(start, end);
  }

  runs(): readonly Run[] | undefined {
    if (!Number.isFinite(this.lo) || !Number.isFinite(this.hi)) return undefined;
    return canonicalRuns(this.forward(this.lo, this.hi));
  }

  private inv?: Linear;
  inverse(): Linear {
    if (this.inv) return this.inv;
    if (this.size !== undefined) this.inv = new RunsLinear(this.runs()!.map(invertRun));
    else if (this.dir === 1) this.inv = new ImplicitLinear(-this.d, 1, this.lo + this.d, this.hi + this.d);
    else this.inv = new ImplicitLinear(this.d, -1, this.d - this.hi + 1, this.d - this.lo + 1);
    return this.inv;
  }
}

// ---- runs ----

/** Runs, canonical; a handful for an offset, a CIGAR or exons. */
export class RunsLinear implements Linear {
  readonly encoding = 'runs';
  readonly list: readonly Run[];
  readonly lo: number;
  readonly hi: number;
  /** The greatest source end of the runs up to each one: monotone, so a binary search finds the first that reaches. */
  private readonly reach: Float64Array;

  constructor(runs: readonly Run[]) {
    this.list = canonicalRuns(runs);
    this.reach = new Float64Array(this.list.length);
    let m = -Infinity;
    this.list.forEach((r, i) => {
      m = Math.max(m, r.from + r.length);
      this.reach[i] = m;
    });
    this.lo = this.list.length ? this.list[0].from : 0;
    this.hi = this.list.length ? m : 0;
  }

  forward(start: number, end: number): Run[] {
    const out: Run[] = [];
    let lo = 0;
    let hi = this.list.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.reach[mid] > start) hi = mid;
      else lo = mid + 1;
    }
    for (let i = lo; i < this.list.length && this.list[i].from < end; i++) {
      const c = clipRun(this.list[i], start, end);
      if (c) out.push(c);
    }
    return out;
  }

  backward(start: number, end: number): Run[] {
    return this.inverse().forward(start, end);
  }

  runs(): readonly Run[] {
    return this.list;
  }

  private inv?: RunsLinear;
  inverse(): Linear {
    if (!this.inv) {
      this.inv = new RunsLinear(this.list.map(invertRun));
      this.inv.inv = this;
    }
    return this.inv;
  }
}

// ---- expanded ----

/** One image a position over `[lo, lo + image.length)`, -1 where none: only for a region asked. */
export class ExpandedLinear implements Linear {
  readonly encoding = 'expanded';
  readonly hi: number;
  constructor(
    readonly lo: number,
    readonly image: Int32Array,
  ) {
    this.hi = lo + image.length;
  }

  /** The runs of `[start, end)`; read one position wider on each side so a lone position keeps its direction. */
  forward(start: number, end: number): Run[] {
    const s = Math.max(start, this.lo);
    const e = Math.min(end, this.hi);
    if (e <= s) return [];
    const wide = this.scan(Math.max(this.lo, s - 1), Math.min(this.hi, e + 1));
    const out: Run[] = [];
    for (const r of wide) {
      const c = clipRun(r, s, e);
      if (c) out.push(c);
    }
    return out;
  }

  private scan(s: number, e: number): Run[] {
    const runs: Run[] = [];
    for (let k = s; k < e; k++) {
      const t = this.image[k - this.lo];
      if (t >= 0) runs.push({ from: k, to: t, length: 1, dir: 1 });
    }
    return canonicalRuns(runs);
  }

  backward(start: number, end: number): Run[] {
    return this.inverse().forward(start, end);
  }

  runs(): readonly Run[] {
    return this.scan(this.lo, this.hi);
  }

  private inv?: Linear;
  inverse(): Linear {
    if (!this.inv) this.inv = new RunsLinear(this.runs().map(invertRun));
    return this.inv;
  }
}

/** A linear map as runs (bounded maps only). */
export function toRuns(l: Linear): RunsLinear {
  const runs = l.runs();
  if (!runs) throw new Error('An unbounded map has no runs');
  return l instanceof RunsLinear ? l : new RunsLinear(runs);
}

/** A linear map expanded over `[lo, hi)`: for tests and the small regions a per-position form is asked for. */
export function expand(l: Linear, lo: number, hi: number): ExpandedLinear {
  const image = new Int32Array(Math.max(0, hi - lo)).fill(-1);
  for (const r of l.forward(lo, hi)) for (let i = 0; i < r.length; i++) image[r.from + i - lo] = r.to + r.dir * i;
  return new ExpandedLinear(lo, image);
}

/** Whether two maps are the same: implicit ones by their arithmetic, others by their canonical runs. */
export function linearsEqual(a: Linear, b: Linear): boolean {
  const ra = a.runs();
  const rb = b.runs();
  if (ra && rb) return runsEqual(ra, rb);
  if (ra || rb) return false;
  const x = a as ImplicitLinear;
  const y = b as ImplicitLinear;
  return x.d === y.d && x.dir === y.dir && x.lo === y.lo && x.hi === y.hi && x.size === y.size;
}

/** Whether a map pairs no position. */
export function isEmptyLinear(l: Linear): boolean {
  const r = l.runs();
  return r ? r.length === 0 : l.hi <= l.lo;
}

// ---- composition ----

/**
 * `b ∘ a`: a's source to b's target. Two implicit maps without wrapping give an implicit one; else runs, each of a's
 * runs carried through b's image of it (an unbounded `a` is composed from the other side, through the inverses).
 */
export function composeLinear(a: Linear, b: Linear): Linear {
  if (a instanceof ImplicitLinear && b instanceof ImplicitLinear && a.size === undefined && b.size === undefined) {
    // b(a(k)) = b.d + b.dir·a.d + b.dir·a.dir·k, for k in a's domain with a(k) in b's
    const [lo, hi] = a.dir === 1 ? [b.lo - a.d, b.hi - a.d] : [a.d - b.hi + 1, a.d - b.lo + 1];
    const l = Math.max(a.lo, lo);
    const h = Math.min(a.hi, hi);
    return new ImplicitLinear(b.d + b.dir * a.d, (a.dir * b.dir) as 1 | -1, l, Math.max(l, h));
  }
  const ra = a.runs();
  if (!ra) {
    if (!b.runs()) throw new Error('Two unbounded maps compose only when implicit');
    return composeLinear(b.inverse(), a.inverse()).inverse();
  }
  const out: Run[] = [];
  for (const r of ra) {
    const [lo, hi] = runTargets(r);
    for (const p of b.forward(lo, hi)) {
      if (r.dir === 1) out.push({ from: r.from + (p.from - r.to), to: p.to, length: p.length, dir: p.dir });
      else {
        const last = p.from + p.length - 1;
        out.push({
          from: r.from + r.to - last,
          to: p.to + p.dir * (p.length - 1),
          length: p.length,
          dir: -p.dir as 1 | -1,
        });
      }
    }
  }
  return new RunsLinear(out);
}

/**
 * A map of residues as one of their nucleotide triplets' places (`3k + j`), for composing with a codon step: a
 * reversed run reverses within the triplet too.
 */
export function scaleLinear(l: Linear, by: 3): Linear {
  if (l instanceof ImplicitLinear) {
    const size = l.size === undefined ? undefined : l.size * by;
    return l.dir === 1
      ? new ImplicitLinear(l.d * by, 1, l.lo * by, l.hi * by, size)
      : new ImplicitLinear(l.d * by + by - 1, -1, l.lo * by, l.hi * by, size);
  }
  return new RunsLinear(
    toRuns(l).list.map((r) => ({
      from: r.from * by,
      to: r.dir === 1 ? r.to * by : r.to * by + by - 1,
      length: r.length * by,
      dir: r.dir,
    })),
  );
}
