/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The approximate matcher of `@molstar/seq-core/algorithms/search` (seqspec:libraries#search): Myers' bit-vector algorithm for
// approximate string matching (JACM 46, 1999), in Hyyrö's blocks of 32 pattern letters (the infix mode of
// `../myers.ts`, which keeps a path and a band for one best alignment; search needs every end, so it has its own loop).
// A forward pass over the text finds the cost of the best alignment ending at each position; a short backward pass from
// each hit finds where it starts. Everything is allocated by the constructor: a search allocates nothing, hits
// included.

import { checkWindow, type HitVisitor, type Matcher, readEnd } from './matchers.js';

/**
 * One text column of Myers' algorithm over the blocks of a pattern: updates the vertical deltas `P` (+1) and `M` (-1)
 * of every block with the match bits `peq[base + block]`, `hin` the horizontal delta entering the top block (0 for an
 * infix search, whose first row is free; +1 for an anchored one). Returns the horizontal delta of the cell at
 * `lastBit` in the last block, the pattern's last row: the change of the score from the previous column.
 */
function column(
  P: Int32Array,
  M: Int32Array,
  peq: Int32Array,
  base: number,
  blocks: number,
  hin: number,
  lastBit: number,
): number {
  let delta = 0;
  for (let b = 0; b < blocks; b++) {
    const pv = P[b];
    const mv = M[b];
    let eq = peq[base + b];
    const xv = eq | mv;
    if (hin < 0) eq |= 1;
    const xh = ((((eq & pv) + pv) | 0) ^ pv) | eq;
    let ph = mv | ~(xh | pv);
    let mh = pv & xh;
    let hout = 0;
    if (b === blocks - 1) delta = (ph & lastBit ? 1 : 0) - (mh & lastBit ? 1 : 0);
    else hout = (ph >>> 31) - (mh >>> 31);
    ph <<= 1;
    mh <<= 1;
    if (hin < 0) mh |= 1;
    else if (hin > 0) ph |= 1;
    P[b] = mh | ~(xv | ph);
    M[b] = ph & xv;
    hin = hout;
  }
  return delta;
}

/** Bit `i` of block `i >>> 5` set in `peq[code * blocks + block]` for each pattern code, and for each code `equal` pairs with it. */
function matchBits(pattern: Uint8Array, equal: Uint8Array | undefined, blocks: number): Int32Array {
  const peq = new Int32Array(256 * blocks);
  for (let i = 0; i < pattern.length; i++) peq[pattern[i] * blocks + (i >>> 5)] |= 1 << (i & 31);
  if (equal) {
    // Pairs match both ways: pattern code a against text code b, and b against a.
    for (let p = 0; p < equal.length; p += 2) {
      for (let i = 0; i < pattern.length; i++) {
        if (pattern[i] === equal[p]) peq[equal[p + 1] * blocks + (i >>> 5)] |= 1 << (i & 31);
        if (pattern[i] === equal[p + 1]) peq[equal[p] * blocks + (i >>> 5)] |= 1 << (i & 31);
      }
    }
  }
  return peq;
}

/**
 * Which ends of the approximate matcher are hits. Each hit is the shortest stretch of its cost ending there.
 * - `"local_best_last"` (the default) and `"local_best_first"`: the local minima of the cost along the text, one hit per
 *   site. A run of ends of equal cost, at most `maxEdits`, with a higher cost on both sides (the start of the text read,
 *   and its end, count as higher) gives one hit, at its last end or at its first; and every end of cost 0 is a hit, so
 *   every exact occurrence is one. The last end reads a mismatch at the pattern's end as a substitution, the first as a
 *   deletion (`RHSWV` in `RHSVVV`: `RHSVV` or `RHSV`).
 * - `"all"`: every end whose cost is at most `maxEdits` (an edit at a site makes several neighboring ends hits).
 */
export type ApproximateEnds = 'local_best_last' | 'local_best_first' | 'all';

export interface ApproximateOptions {
  /** Pairs of codes `[a0, b0, a1, b1, …]` that also match, both ways, besides equal codes (N against every nucleotide). */
  readonly equal?: Uint8Array;
  /** Which ends are hits; default `"local_best_last"`. */
  readonly ends?: ApproximateEnds;
}

class ApproximateMatcher implements Matcher {
  readonly length: number;
  readonly m: number;
  readonly k: number;
  readonly all: boolean;
  readonly last: boolean;
  readonly blocks: number;
  readonly lastBit: number;
  readonly peq: Int32Array;
  readonly peqBack: Int32Array;
  readonly P: Int32Array;
  readonly M: Int32Array;
  readonly backP: Int32Array;
  readonly backM: Int32Array;

  constructor(pattern: Uint8Array, maxEdits: number, options: ApproximateOptions | undefined) {
    const m = pattern.length;
    const equal = options?.equal;
    const ends = options?.ends ?? 'local_best_last';
    if (!Number.isInteger(maxEdits) || maxEdits < 0) {
      throw new RangeError(`approximateMatcher: maxEdits must be an integer of at least 0; got ${maxEdits}`);
    }
    if (m > 0 && maxEdits >= m) {
      throw new RangeError(
        `approximateMatcher: ${maxEdits} edits in a pattern of ${m} letters would match anything; use at most ${m - 1}`,
      );
    }
    if (equal && equal.length % 2 !== 0)
      throw new RangeError('approximateMatcher: equal is pairs of codes, an even count');
    if (ends !== 'local_best_last' && ends !== 'local_best_first' && ends !== 'all') {
      throw new RangeError(
        `approximateMatcher: ends is "local_best_last", "local_best_first" or "all"; got ${String(ends)}`,
      );
    }
    this.m = m;
    this.k = maxEdits;
    this.all = ends === 'all';
    this.last = ends === 'local_best_last';
    this.length = m === 0 ? 0 : m + maxEdits;
    this.blocks = (m + 31) >>> 5;
    this.lastBit = m === 0 ? 0 : 1 << ((m - 1) & 31);
    this.peq = matchBits(pattern, equal, this.blocks);
    this.peqBack = matchBits(pattern.slice().reverse(), equal, this.blocks);
    this.P = new Int32Array(this.blocks);
    this.M = new Int32Array(this.blocks);
    this.backP = new Int32Array(this.blocks);
    this.backM = new Int32Array(this.blocks);
  }

  find(seq: Uint8Array, from: number, to: number, cyclic: boolean, visit: HitVisitor): number {
    checkWindow(seq.length, from, to);
    if (this.m === 0 || from >= to) return 0;
    const { m, k, all, last, blocks, lastBit, peq, P, M } = this;
    const n = seq.length;
    // One past the last end of a hit that starts in the window. The text read runs to `cap` (a lap past `to` on a
    // circle): the local-best policies read on past `limit` while a run of equal cost is pending, to see whether it is
    // a minimum (and, "last", where it ends).
    const limit = readEnd(this.length, n, to, cyclic);
    const cap = cyclic ? to + n : to;
    P.fill(-1);
    M.fill(0);
    let score = m;
    let before = m;
    let pendingEnd = -1;
    let pendingScore = 0;
    let count = 0;
    for (let u = from; u < limit || (pendingEnd >= 0 && u < cap); u++) {
      score += column(P, M, peq, seq[u >= n ? u - n : u] * blocks, blocks, 0, lastBit);
      if (all) {
        if (score <= k) count += this.report(seq, from, to, u, score, visit);
        continue;
      }
      // A pending run ends: a minimum when the cost rises, nothing when it falls; while it lasts, "last" follows it.
      if (pendingEnd >= 0) {
        if (score === pendingScore) {
          if (last) pendingEnd = u;
        } else {
          if (score > pendingScore) count += this.report(seq, from, to, pendingEnd, pendingScore, visit);
          pendingEnd = -1;
        }
      }
      if (score === 0) count += this.report(seq, from, to, u, 0, visit);
      else if (score <= k && score < before && u < limit) {
        pendingEnd = u;
        pendingScore = score;
      }
      before = score;
    }
    if (pendingEnd >= 0) count += this.report(seq, from, to, pendingEnd, pendingScore, visit);
    return count;
  }

  /**
   * The hit ending at `end` (inclusive, in the unrolled text) with cost `cost`: the shortest stretch ending there that
   * costs `cost`, found by running the reversed pattern leftwards from `end` with the text's start anchored. Reports it
   * when it starts in the window; returns 1 then, else 0.
   */
  private report(seq: Uint8Array, from: number, to: number, end: number, cost: number, visit: HitVisitor): number {
    const { m, blocks, lastBit, peqBack, backP, backM } = this;
    const n = seq.length;
    backP.fill(-1);
    backM.fill(0);
    let score = m;
    // At most a lap, and the stretches of an alignment of cost `cost`.
    const longest = Math.min(m + cost, end + 1 - from, n);
    for (let t = 1; t <= longest; t++) {
      const u = end - t + 1;
      score += column(backP, backM, peqBack, seq[u >= n ? u - n : u] * blocks, blocks, 1, lastBit);
      if (score === cost) {
        if (u >= to) return 0;
        visit(u, end + 1 > n ? end + 1 - n : end + 1, cost);
        return 1;
      }
    }
    return 0;
  }
}

/**
 * Hits where `pattern` occurs with at most `maxEdits` edits (substitutions, insertions, deletions), in the infix mode
 * of Myers' bit-vector algorithm. For each end position the cost is the cheapest alignment of the pattern to a stretch
 * ending there; `options.ends` says which ends are hits (`ApproximateEnds`: by default the local minima of the cost,
 * one hit per site at the last end of a run of equal cost, every exact occurrence included). A hit is the shortest stretch of its cost ending there; `score`
 * is the cost. `options.equal`: code pairs that also match. `maxEdits` is at most the pattern's length minus one (more
 * would match anything). With a pattern that fits a word (32 letters) a text letter costs a handful of integer
 * operations; each further 32 letters add one pass; a hit adds a pass over at most `length` letters.
 */
export function approximateMatcher(pattern: Uint8Array, maxEdits: number, options?: ApproximateOptions): Matcher {
  return new ApproximateMatcher(pattern, maxEdits, options);
}
