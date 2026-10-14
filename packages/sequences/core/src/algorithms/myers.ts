/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `myers` aligner (A2; seqspec:format#mappings): unit-cost edit distance with an alignment path, by Myers'
// bit-vector algorithm as edlib implements it, with Ukkonen's band, a traceback for small problems and Hirschberg's
// split for large ones.
//
// Ported from edlib, https://github.com/Martinsos/edlib, commit 0ddc23ea06abd95ab6227442acb4c27bd7607b68,
// edlib/src/edlib.cpp (edlibAlign, buildPeq, calculateBlock, myersCalcEditDistanceSemiGlobal,
// myersCalcEditDistanceNW, obtainAlignment, obtainAlignmentTraceback, obtainAlignmentHirschberg, transformSequences).
// What differs from edlib: VENDORED.md, "A2". edlib's license:
//
//   The MIT License (MIT)
//
//   Copyright (c) 2014 Martin Šošić
//
//   Permission is hereby granted, free of charge, to any person obtaining a copy of
//   this software and associated documentation files (the "Software"), to deal in
//   the Software without restriction, including without limitation the rights to
//   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//   the Software, and to permit persons to whom the Software is furnished to do so,
//   subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be included in all
//   copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//   COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//   IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// Words are 32 bits (int32 in an Int32Array, bitwise operators, `>>>`) where edlib's are 64 (`uint64_t`): a block
// covers 32 query letters. A word's high bit is set exactly when it is negative, so `x < 0` tests edlib's
// `x & HIGH_BIT_MASK`. edlib's query is `a` and its target `b`; its INSERT (a letter only in the query) is our `I`
// and its DELETE (a letter only in the target) our `D` (seqspec:model#coordinate-maps), so its path reads as our CIGAR
// unchanged.

import { type Aligner, type AlignerParams, alignable, readOptions } from './aligner.js';

const WORD_SIZE = 32;
/** edlib's EDLIB_EDOP_*: INSERT a letter only in `a` (our I), DELETE a letter only in `b` (our D). */
const EDOP = { MATCH: 0, INSERT: 1, DELETE: 2, MISMATCH: 3 } as const;
/** edlib's modes: NW global, SHW prefix (gaps after `a` free), HW infix (gaps before and after `a` free). */
const MODE = { NW: 0, SHW: 1, HW: 2 } as const;
/** Every STRONG_REDUCE_NUM columns the band is reduced in the more expensive, cell-by-cell way. */
const STRONG_REDUCE_NUM = 2048;
/** Below this many bytes of stored columns the path comes from a traceback; above it, Hirschberg's split (edlib's 1 MB). */
const TRACEBACK_BYTES = 1024 * 1024;
/** The most words (and columns) a traceback can store under TRACEBACK_BYTES: 12 bytes per block, 8 per column. */
const TRACEBACK_WORDS = Math.ceil(TRACEBACK_BYTES / 12);
const TRACEBACK_COLUMNS = TRACEBACK_BYTES / 8;
/** `nw`'s result when it stopped at `targetStopPosition` and stored that column. */
const STOPPED = -2;

/** Buffers for one alignment, allocated once and sized for the whole problem; sub-problems use their prefixes. */
interface Work {
  readonly alphabetLength: number;
  /** alphabetLength²: 1 where two letters match (themselves, and the additional equalities). */
  readonly equal: Uint8Array;
  readonly peq: Int32Array;
  readonly rPeq: Int32Array;
  /** The blocks of the current column (edlib's `Block`: Pv, Mv and the score of the block's last cell). */
  readonly P: Int32Array;
  readonly M: Int32Array;
  readonly score: Int32Array;
  /** edlib's AlignmentData: stored columns for the traceback, or the one column `nw` stops at. */
  readonly dataP: Int32Array;
  readonly dataM: Int32Array;
  readonly dataScore: Int32Array;
  readonly firstBlocks: Int32Array;
  readonly lastBlocks: Int32Array;
  /** Hirschberg's middle columns, unwrapped into scores. */
  readonly left: Int32Array;
  readonly right: Int32Array;
  /** The path, one EDOP per column, appended at `n`. */
  readonly ops: Uint8Array;
  n: number;
  /** `semiGlobal`'s first and last end position in the target with the best score. */
  first: number;
  last: number;
}

/** An alignment found by `myersAlign`. */
export interface MyersResult {
  /** The edit distance: substitutions, insertions and deletions, each costing 1 (free end gaps in `b` excluded). */
  readonly distance: number;
  /** From `a` to `b`, covering both: the unmatched ends of `b` in `prefix` and `infix` are `D`. */
  readonly cigar: string;
  /** Where in `b` the aligned part starts and ends (inclusive); `end` is -1 when `a` aligns before `b` entirely. */
  readonly start: number;
  readonly end: number;
}

export type MyersMode = 'global' | 'prefix' | 'infix';

/**
 * The unit-cost alignment of `a` to `b`: `global` end to end, `prefix` `a` against the start of `b`, `infix` `a`
 * anywhere inside `b`. Null when the distance exceeds `maxDistance` (null: no bound). `equalities` are letter pairs
 * that match besides equal letters, e.g. `["NA", "NC"]`. Letters are case-insensitive; `a` and `b` are letters only and
 * not empty. Throws a `RangeError` on any other input.
 */
export function myersAlign(
  a: string,
  b: string,
  mode: MyersMode = 'global',
  maxDistance: number | null = null,
  equalities: readonly string[] = [],
): MyersResult | null {
  const queryLength = a.length;
  const targetLength = b.length;
  if (queryLength === 0 || targetLength === 0) throw new RangeError('myers: empty sequence');
  // Letters only: `& 127` below would alias a non-ASCII character to an ASCII one.
  if (!alignable(a, b)) throw new RangeError('myers: sequences must be letters A-Z only');
  if (maxDistance !== null && !(Number.isInteger(maxDistance) && maxDistance >= 0)) {
    throw new RangeError(`myers: maxDistance must be null or an integer of at least 0; got ${maxDistance}`);
  }
  for (const pair of equalities) {
    if (!/^[A-Za-z]{2}$/.test(pair)) throw new RangeError(`myers: an equality is two letters; got "${pair}"`);
  }
  const m = mode === 'prefix' ? MODE.SHW : mode === 'infix' ? MODE.HW : MODE.NW;

  // transformSequences: letters to their index in the alphabet, in order of appearance.
  const index = new Int16Array(128).fill(-1);
  let alphabetLength = 0;
  const encode = (s: string) => {
    const out = new Uint8Array(s.length);
    for (let i = 0; i < s.length; i++) {
      const c = s.charCodeAt(i) & 127;
      if (index[c] < 0) index[c] = alphabetLength++;
      out[i] = index[c];
    }
    return out;
  };
  const query = encode(a.toUpperCase());
  const target = encode(b.toUpperCase());

  const maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const W = maxNumBlocks * WORD_SIZE - queryLength;
  const data = Math.max(maxNumBlocks, Math.min(maxNumBlocks * targetLength, TRACEBACK_WORDS));
  const w: Work = {
    alphabetLength,
    equal: new Uint8Array(alphabetLength * alphabetLength),
    peq: new Int32Array(alphabetLength * maxNumBlocks),
    rPeq: new Int32Array(alphabetLength * maxNumBlocks),
    P: new Int32Array(maxNumBlocks),
    M: new Int32Array(maxNumBlocks),
    score: new Int32Array(maxNumBlocks),
    dataP: new Int32Array(data),
    dataM: new Int32Array(data),
    dataScore: new Int32Array(data),
    firstBlocks: new Int32Array(Math.min(targetLength, TRACEBACK_COLUMNS)),
    lastBlocks: new Int32Array(Math.min(targetLength, TRACEBACK_COLUMNS)),
    left: new Int32Array(maxNumBlocks * WORD_SIZE),
    right: new Int32Array(maxNumBlocks * WORD_SIZE),
    ops: new Uint8Array(queryLength + targetLength),
    n: 0,
    first: -1,
    last: -1,
  };
  // EqualityDefinition: each letter equals itself; a pair naming a letter neither sequence has is ignored.
  for (let i = 0; i < alphabetLength; i++) w.equal[i * alphabetLength + i] = 1;
  for (const pair of equalities) {
    const x = index[pair.toUpperCase().charCodeAt(0) & 127];
    const y = index[pair.toUpperCase().charCodeAt(1) & 127];
    if (x >= 0 && y >= 0) w.equal[x * alphabetLength + y] = w.equal[y * alphabetLength + x] = 1;
  }
  buildPeq(w, w.peq, query, queryLength);

  // Without a bound, k starts at a word and doubles until a solution is found.
  const dynamicK = maxDistance === null;
  let k = maxDistance ?? WORD_SIZE;
  let distance: number;
  do {
    distance =
      m === MODE.NW
        ? nw(w, w.peq, W, maxNumBlocks, queryLength, target, targetLength, k, false, -1)
        : semiGlobal(w, w.peq, W, maxNumBlocks, queryLength, target, targetLength, k, m);
    k *= 2;
  } while (dynamicK && distance < 0);
  if (distance < 0) return null;

  // Start and end locations, the alignment for the first optimal end location only (as edlib).
  const end = m === MODE.NW ? targetLength - 1 : w.first;
  let start = 0;
  const rQuery = reverse(query);
  const rTarget = reverse(target);
  if (m === MODE.HW && end !== -1) {
    // The start from the reversed query against the reversed target ending at `end`; its last best location, so
    // the path does not start with insertions when it can start with mismatches.
    buildPeq(w, w.rPeq, rQuery, queryLength);
    semiGlobal(
      w,
      w.rPeq,
      W,
      maxNumBlocks,
      queryLength,
      rTarget.subarray(targetLength - end - 1),
      end + 1,
      distance,
      MODE.SHW,
    );
    start = end - w.last;
  }
  w.ops.fill(EDOP.DELETE, 0, start);
  w.n = start;
  const alnTargetLength = end - start + 1;
  obtainAlignment(
    w,
    query,
    rQuery,
    queryLength,
    target.subarray(start),
    rTarget.subarray(targetLength - 1 - end),
    alnTargetLength,
    distance,
  );
  w.ops.fill(EDOP.DELETE, w.n, w.n + targetLength - 1 - end);
  w.n += targetLength - 1 - end;
  return { distance, cigar: cigarOfOps(w.ops, w.n), start, end };
}

const params: AlignerParams = {
  mode: {
    type: 'choice',
    choices: ['global', 'prefix', 'infix'],
    nullable: false,
    default: 'global',
    description:
      'global: end to end; prefix: `a` against the start of `b`; infix: `a` found anywhere inside `b` (a construct ' +
      'in its UniProt sequence). The unmatched ends of `b` become D.',
  },
  max_distance: {
    type: 'int',
    nullable: true,
    min: 0,
    default: null,
    description: 'The largest edit distance accepted; beyond it the sequences are not aligned. null: no bound.',
  },
  equalities: {
    type: 'str',
    pattern: '[A-Za-z]{2}([ ,]+[A-Za-z]{2})*',
    nullable: true,
    default: null,
    description: 'Letter pairs that match besides equal letters, separated by spaces, e.g. "NA NC NG NT". null: none.',
  },
};

export const myers: Aligner = {
  name: 'myers',
  description:
    "Unit-cost edit distance (Myers' bit-vector algorithm, ported from edlib), for near-identical sequences: " +
    "fast and exact, but gaps cost as much as substitutions. Memory: a bit per cell up to 1 MB, then Hirschberg's split.",
  params,
  async align(a, b, options) {
    const o = readOptions('myers', params, options);
    if (!alignable(a, b)) return null;
    const pairs = o.equalities === null ? [] : (o.equalities as string).split(/[ ,]+/);
    return myersAlign(a, b, o.mode as MyersMode, o.max_distance as number | null, pairs)?.cigar ?? null;
  },
};

/** A CIGAR from edit operations, matches and mismatches both as `M`. */
function cigarOfOps(ops: Uint8Array, n: number): string {
  let cigar = '';
  for (let k = 0; k < n; ) {
    const op = ops[k] === EDOP.INSERT ? 'I' : ops[k] === EDOP.DELETE ? 'D' : 'M';
    let r = k + 1;
    while (r < n && (ops[r] === EDOP.INSERT ? 'I' : ops[r] === EDOP.DELETE ? 'D' : 'M') === op) r++;
    cigar += `${r - k}${op}`;
    k = r;
  }
  return cigar;
}

function ceilDiv(x: number, y: number): number {
  return x % y ? ((x / y) | 0) + 1 : (x / y) | 0;
}

function reverse(seq: Uint8Array): Uint8Array {
  const out = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) out[i] = seq[seq.length - 1 - i];
  return out;
}

/**
 * buildPeq: bit i of `out[s * maxNumBlocks + b]` is set when letter i of the query's block b equals symbol s. The
 * query is padded at its end with letters that match anything. (edlib's extra wildcard row is never read: dropped.)
 */
function buildPeq(w: Work, out: Int32Array, query: Uint8Array, queryLength: number): void {
  const maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const { alphabetLength, equal } = w;
  for (let symbol = 0; symbol < alphabetLength; symbol++) {
    for (let b = 0; b < maxNumBlocks; b++) {
      let v = 0;
      for (let r = (b + 1) * WORD_SIZE - 1; r >= b * WORD_SIZE; r--) {
        v = (v << 1) | (r >= queryLength || equal[query[r] * alphabetLength + symbol] ? 1 : 0);
      }
      out[symbol * maxNumBlocks + b] = v;
    }
  }
}

/**
 * calculateBlock (Myers' Advance_Block): advances block `i` of `P`/`M` by one column whose match bits are `eq`, `hin`
 * the horizontal delta entering it at the top (+1, 0 or -1); writes the new Pv/Mv in place, returns the delta out.
 */
function advance(P: Int32Array, M: Int32Array, i: number, eq: number, hin: number): number {
  const pv = P[i];
  const mv = M[i];
  const hinIsNeg = (hin >> 2) & 1;
  const xv = eq | mv;
  eq |= hinIsNeg;
  const xh = (((eq & pv) + pv) ^ pv) | eq;
  let ph = mv | ~(xh | pv);
  let mh = pv & xh;
  const hout = (ph >>> 31) - (mh >>> 31);
  ph = (ph << 1) | ((hin + 1) >> 1);
  mh = (mh << 1) | hinIsNeg;
  P[i] = mh | ~(xv | ph);
  M[i] = ph & xv;
  return hout;
}

/** The score of cell `i` of a block counted from its bottom (cell 0 has `score`): getBlockCellValues(block)[i]. */
function cellFromBottom(p: number, m: number, score: number, i: number): number {
  for (let j = 0; j < i; j++, p <<= 1, m <<= 1) score += (m >>> 31) - (p >>> 31);
  return score;
}

/** allBlockCellsLarger: whether every cell of the block scores more than k. */
function allLarger(p: number, m: number, score: number, k: number): boolean {
  for (let j = 0; j < WORD_SIZE; j++, p <<= 1, m <<= 1) {
    if (score <= k) return false;
    score += (m >>> 31) - (p >>> 31);
  }
  return true;
}

/**
 * The strong reduction's test for one block of `nw`: whether a cell of its `numCells` real cells (the bottom one in
 * row `bottom`) scores at most k and lies inside the band: above its lower edge (`lower`), or below its upper edge.
 */
function anyCellInBand(
  p: number,
  m: number,
  score: number,
  numCells: number,
  bottom: number,
  k: number,
  offset: number,
  lower: boolean,
): boolean {
  let r = bottom;
  for (let i = 0; i < WORD_SIZE; i++, p <<= 1, m <<= 1) {
    if (i >= WORD_SIZE - numCells) {
      if (score <= k && (lower ? r <= k - score + offset + 1 : r >= score - k + offset)) return true;
      r--;
    }
    score += (m >>> 31) - (p >>> 31);
  }
  return false;
}

/**
 * myersCalcEditDistanceSemiGlobal (modes SHW and HW): the best score at most k over every end position in the target,
 * -1 when none; `w.first` and `w.last` are the first and last end positions with that score.
 */
function semiGlobal(
  w: Work,
  peq: Int32Array,
  W: number,
  maxNumBlocks: number,
  queryLength: number,
  target: Uint8Array,
  targetLength: number,
  k: number,
  mode: number,
): number {
  const { P, M, score } = w;
  w.first = w.last = -1;
  // The Ukkonen band: its first and last block (edlib's `bl` always equals lastBlock).
  let firstBlock = 0;
  let lastBlock = Math.min(ceilDiv(k + 1, WORD_SIZE), maxNumBlocks) - 1;
  // For HW, the solution is never larger than the query's length.
  if (mode === MODE.HW) k = Math.min(queryLength, k);
  for (let b = 0; b <= lastBlock; b++) {
    score[b] = (b + 1) * WORD_SIZE;
    P[b] = -1;
    M[b] = 0;
  }
  let bestScore = -1;
  // HW: a gap before the query is free (hin 0 at the top).
  const startHout = mode === MODE.HW ? 0 : 1;
  for (let c = 0; c < targetLength; c++) {
    const pc = target[c] * maxNumBlocks;
    let hout = startHout;
    for (let b = firstBlock; b <= lastBlock; b++) {
      hout = advance(P, M, b, peq[pc + b], hout);
      score[b] += hout;
    }

    // Adjust the band (Ukkonen): one more block when the last one's score is not too big.
    if (
      lastBlock < maxNumBlocks - 1 &&
      score[lastBlock] - hout <= k &&
      ((peq[pc + lastBlock + 1] & 1) !== 0 || hout < 0)
    ) {
      lastBlock++;
      P[lastBlock] = -1;
      M[lastBlock] = 0;
      score[lastBlock] = score[lastBlock - 1] - hout + WORD_SIZE + advance(P, M, lastBlock, peq[pc + lastBlock], hout);
    } else {
      while (lastBlock >= firstBlock && score[lastBlock] >= k + WORD_SIZE) lastBlock--;
    }
    if (c % STRONG_REDUCE_NUM === 0) {
      while (lastBlock >= 0 && lastBlock >= firstBlock && allLarger(P[lastBlock], M[lastBlock], score[lastBlock], k)) {
        lastBlock--;
      }
    }
    // For HW the first block is always a candidate (the top boundary is 0): the band never ends before the last column.
    if (mode === MODE.HW && lastBlock === -1) lastBlock++;
    if (mode !== MODE.HW) {
      while (firstBlock <= lastBlock && score[firstBlock] >= k + WORD_SIZE) firstBlock++;
      if (c % STRONG_REDUCE_NUM === 0) {
        while (firstBlock <= lastBlock && allLarger(P[firstBlock], M[firstBlock], score[firstBlock], k)) firstBlock++;
      }
    }
    if (lastBlock < firstBlock) return bestScore;

    // The best score: the score in column c is the one of column c - W (the padding matches diagonally).
    if (lastBlock === maxNumBlocks - 1) {
      const colScore = score[lastBlock];
      if (colScore <= k && (bestScore === -1 || colScore <= bestScore)) {
        if (colScore !== bestScore) {
          k = bestScore = colScore;
          w.first = c - W;
        }
        w.last = c - W;
      }
    }
  }
  // The last W columns' scores, from the last column's padded cells.
  if (lastBlock === maxNumBlocks - 1) {
    for (let i = 0; i < W; i++) {
      const colScore = cellFromBottom(P[lastBlock], M[lastBlock], score[lastBlock], i + 1);
      if (colScore <= k && (bestScore === -1 || colScore <= bestScore)) {
        if (colScore !== bestScore) {
          k = bestScore = colScore;
          w.first = targetLength - W + i;
        }
        w.last = targetLength - W + i;
      }
    }
  }
  return bestScore;
}

/**
 * myersCalcEditDistanceNW: the global edit distance when at most k, -1 otherwise. `findAlignment` stores every
 * column for the traceback; `targetStopPosition` > -1 stops at that column, stores it alone and returns STOPPED.
 */
function nw(
  w: Work,
  peq: Int32Array,
  W: number,
  maxNumBlocks: number,
  queryLength: number,
  target: Uint8Array,
  targetLength: number,
  k: number,
  findAlignment: boolean,
  targetStopPosition: number,
): number {
  if (k < Math.abs(targetLength - queryLength)) return -1;
  const { P, M, score } = w;
  k = Math.min(k, Math.max(queryLength, targetLength));
  let firstBlock = 0;
  let lastBlock =
    Math.min(maxNumBlocks, ceilDiv(Math.min(k, ((k + queryLength - targetLength) / 2) | 0) + 1, WORD_SIZE)) - 1;
  for (let b = 0; b <= lastBlock; b++) {
    score[b] = (b + 1) * WORD_SIZE;
    P[b] = -1;
    M[b] = 0;
  }
  // The band's edges below are in rows; `diag` is their shared offset for column c.
  for (let c = 0; c < targetLength; c++) {
    const pc = target[c] * maxNumBlocks;
    let hout = 1;
    for (let b = firstBlock; b <= lastBlock; b++) {
      hout = advance(P, M, b, peq[pc + b], hout);
      score[b] += hout;
    }
    // Update k at the end of the column; in the last block the score is from W cells left and up.
    k = Math.min(
      k,
      score[lastBlock] +
        Math.max(targetLength - c - 1, queryLength - ((1 + lastBlock) * WORD_SIZE - 1) - 1) +
        (lastBlock === maxNumBlocks - 1 ? W : 0),
    );
    const diag = c + queryLength - targetLength;

    // The last block: calculate the next one when it is not beneath the band, then drop blocks out of the band.
    if (
      lastBlock + 1 < maxNumBlocks &&
      !((lastBlock + 1) * WORD_SIZE - 1 > k - score[lastBlock] + 2 * WORD_SIZE - 2 + diag)
    ) {
      lastBlock++;
      P[lastBlock] = -1;
      M[lastBlock] = 0;
      const newHout = advance(P, M, lastBlock, peq[pc + lastBlock], hout);
      score[lastBlock] = score[lastBlock - 1] - hout + WORD_SIZE + newHout;
      hout = newHout;
    }
    // (edlib: "Does not work if do not put +1! Why???")
    while (
      lastBlock >= firstBlock &&
      (score[lastBlock] >= k + WORD_SIZE ||
        (lastBlock + 1) * WORD_SIZE - 1 > k - score[lastBlock] + 2 * WORD_SIZE - 2 + diag + 1)
    ) {
      lastBlock--;
    }
    // The first block: advance while outside the band.
    while (
      firstBlock <= lastBlock &&
      (score[firstBlock] >= k + WORD_SIZE || (firstBlock + 1) * WORD_SIZE - 1 < score[firstBlock] - k + diag)
    ) {
      firstBlock++;
    }
    if (c % STRONG_REDUCE_NUM === 0) {
      while (lastBlock >= firstBlock) {
        const numCells = lastBlock === maxNumBlocks - 1 ? WORD_SIZE - W : WORD_SIZE;
        const bottom = lastBlock * WORD_SIZE + numCells - 1;
        if (anyCellInBand(P[lastBlock], M[lastBlock], score[lastBlock], numCells, bottom, k, diag, true)) break;
        lastBlock--;
      }
      while (firstBlock <= lastBlock) {
        const numCells = firstBlock === maxNumBlocks - 1 ? WORD_SIZE - W : WORD_SIZE;
        const bottom = firstBlock * WORD_SIZE + numCells - 1;
        if (anyCellInBand(P[firstBlock], M[firstBlock], score[firstBlock], numCells, bottom, k, diag, false)) break;
        firstBlock++;
      }
    }
    if (lastBlock < firstBlock) return -1;

    if (findAlignment) {
      const at = maxNumBlocks * c;
      for (let b = firstBlock; b <= lastBlock; b++) {
        w.dataP[at + b] = P[b];
        w.dataM[at + b] = M[b];
        w.dataScore[at + b] = score[b];
      }
      w.firstBlocks[c] = firstBlock;
      w.lastBlocks[c] = lastBlock;
    }
    if (c === targetStopPosition) {
      for (let b = firstBlock; b <= lastBlock; b++) {
        w.dataP[b] = P[b];
        w.dataM[b] = M[b];
        w.dataScore[b] = score[b];
      }
      w.firstBlocks[0] = firstBlock;
      w.lastBlocks[0] = lastBlock;
      return STOPPED;
    }
  }
  if (lastBlock === maxNumBlocks - 1) {
    // The best score is W cells above the last block's bottom: the query is padded with W cells.
    const bestScore = cellFromBottom(P[lastBlock], M[lastBlock], score[lastBlock], W);
    if (bestScore <= k) return bestScore;
  }
  return -1;
}

/** obtainAlignment: appends an optimal path (one of score `bestScore`) to `w.ops`. */
function obtainAlignment(
  w: Work,
  query: Uint8Array,
  rQuery: Uint8Array,
  queryLength: number,
  target: Uint8Array,
  rTarget: Uint8Array,
  targetLength: number,
  bestScore: number,
): void {
  if (queryLength === 0 || targetLength === 0) {
    w.ops.fill(queryLength === 0 ? EDOP.DELETE : EDOP.INSERT, w.n, w.n + queryLength + targetLength);
    w.n += queryLength + targetLength;
    return;
  }
  const maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const W = maxNumBlocks * WORD_SIZE - queryLength;
  // The traceback when its stored columns (two words and a score per block, two ints per column) stay under 1 MB.
  if ((2 * 4 + 4) * maxNumBlocks * targetLength + 2 * 4 * targetLength < TRACEBACK_BYTES) {
    buildPeq(w, w.peq, query, queryLength);
    nw(w, w.peq, W, maxNumBlocks, queryLength, target, targetLength, bestScore, true, -1);
    traceback(w, queryLength, targetLength, bestScore);
  } else {
    hirschberg(w, query, rQuery, queryLength, target, rTarget, targetLength, bestScore);
  }
}

/**
 * obtainAlignmentTraceback: one optimal path, back through the columns `nw` stored, from the bottom right cell.
 * Preference at each cell: up (I), left (D), then diagonal (M).
 */
function traceback(w: Work, queryLength: number, targetLength: number, bestScore: number): void {
  const maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const W = maxNumBlocks * WORD_SIZE - queryLength;
  const { dataP: Ps, dataM: Ms, dataScore: scores, firstBlocks, lastBlocks, ops } = w;
  const inBand = (c: number, b: number) => c > 0 && b >= firstBlocks[c - 1] && b <= lastBlocks[c - 1];
  const begin = w.n;
  let n = w.n;
  let c = targetLength - 1; // the column
  let b = maxNumBlocks - 1; // the block in the column
  let currScore = bestScore;
  let lScore = -1; // left
  let uScore = -1; // upper
  let ulScore = -1; // upper left
  let currP = Ps[c * maxNumBlocks + b];
  let currM = Ms[c * maxNumBlocks + b];
  let thereIsLeftBlock = inBand(c, b);
  let lP = 0;
  let lM = 0;
  if (thereIsLeftBlock) {
    lP = Ps[(c - 1) * maxNumBlocks + b];
    lM = Ms[(c - 1) * maxNumBlocks + b];
  }
  currP <<= W;
  currM <<= W;
  let blockPos = WORD_SIZE - W - 1; // the current cell in the block
  // The new left block after moving left or diagonally; at column 0 it is the boundary.
  const setLeft = () => {
    if (inBand(c, b)) {
      thereIsLeftBlock = true;
      lP = Ps[(c - 1) * maxNumBlocks + b];
      lM = Ms[(c - 1) * maxNumBlocks + b];
    } else if (c === 0) {
      thereIsLeftBlock = true;
      lScore = b * WORD_SIZE + blockPos + 1;
      ulScore = lScore - 1;
    } else {
      thereIsLeftBlock = false;
    }
  };

  while (true) {
    if (c === 0) {
      thereIsLeftBlock = true;
      lScore = b * WORD_SIZE + blockPos + 1;
      ulScore = lScore - 1;
    }
    // Scores of the left, upper left and upper cells.
    if (lScore === -1 && thereIsLeftBlock) {
      lScore = scores[(c - 1) * maxNumBlocks + b];
      for (let i = 0; i < WORD_SIZE - blockPos - 1; i++, lP <<= 1, lM <<= 1) lScore += (lM >>> 31) - (lP >>> 31);
    }
    if (ulScore === -1) {
      if (lScore !== -1) {
        ulScore = lScore + (lM >>> 31) - (lP >>> 31);
      } else if (inBand(c, b - 1)) {
        // The upper left cell is the last cell of a block whose left neighbor is not in the band.
        ulScore = scores[(c - 1) * maxNumBlocks + b - 1];
      }
    }
    if (uScore === -1) {
      uScore = currScore + (currM >>> 31) - (currP >>> 31);
      currP <<= 1;
      currM <<= 1;
    }

    if (uScore !== -1 && uScore + 1 === currScore) {
      // Up: a letter only in the query.
      currScore = uScore;
      lScore = ulScore;
      uScore = ulScore = -1;
      if (blockPos === 0) {
        if (b === 0) {
          // Only boundary cells above: left to the end.
          ops[n++] = EDOP.INSERT;
          for (let i = 0; i < c + 1; i++) ops[n++] = EDOP.DELETE;
          break;
        }
        blockPos = WORD_SIZE - 1;
        b--;
        currP = Ps[c * maxNumBlocks + b];
        currM = Ms[c * maxNumBlocks + b];
        thereIsLeftBlock = inBand(c, b);
        if (thereIsLeftBlock) {
          lP = Ps[(c - 1) * maxNumBlocks + b];
          lM = Ms[(c - 1) * maxNumBlocks + b];
        }
      } else {
        blockPos--;
        lP <<= 1;
        lM <<= 1;
      }
      ops[n++] = EDOP.INSERT;
    } else if (lScore !== -1 && lScore + 1 === currScore) {
      // Left: a letter only in the target.
      currScore = lScore;
      uScore = ulScore;
      lScore = ulScore = -1;
      c--;
      if (c === -1) {
        // Only boundary cells to the left: up to the end.
        ops[n++] = EDOP.DELETE;
        const numUp = b * WORD_SIZE + blockPos + 1;
        for (let i = 0; i < numUp; i++) ops[n++] = EDOP.INSERT;
        break;
      }
      currP = lP;
      currM = lM;
      setLeft();
      ops[n++] = EDOP.DELETE;
    } else if (ulScore !== -1) {
      // Diagonal: a match or a mismatch.
      const moveCode = ulScore === currScore ? EDOP.MATCH : EDOP.MISMATCH;
      currScore = ulScore;
      uScore = lScore = ulScore = -1;
      c--;
      if (c === -1) {
        ops[n++] = moveCode;
        const numUp = b * WORD_SIZE + blockPos;
        for (let i = 0; i < numUp; i++) ops[n++] = EDOP.INSERT;
        break;
      }
      if (blockPos === 0) {
        if (b === 0) {
          ops[n++] = moveCode;
          for (let i = 0; i < c + 1; i++) ops[n++] = EDOP.DELETE;
          break;
        }
        blockPos = WORD_SIZE - 1;
        b--;
        currP = Ps[c * maxNumBlocks + b];
        currM = Ms[c * maxNumBlocks + b];
      } else {
        blockPos--;
        currP = lP << 1;
        currM = lM << 1;
      }
      setLeft();
      ops[n++] = moveCode;
    } else {
      break;
    }
  }
  ops.subarray(begin, n).reverse();
  w.n = n;
}

/**
 * obtainAlignmentHirschberg: splits the target in half, finds the row where an optimal path crosses from the left
 * half's last column to the right half's first (from the forward and the reversed problem), and solves both quarters.
 */
function hirschberg(
  w: Work,
  query: Uint8Array,
  rQuery: Uint8Array,
  queryLength: number,
  target: Uint8Array,
  rTarget: Uint8Array,
  targetLength: number,
  bestScore: number,
): void {
  const maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const W = maxNumBlocks * WORD_SIZE - queryLength;
  const leftHalfWidth = targetLength >> 1;
  const rightHalfWidth = targetLength - leftHalfWidth;
  const { dataP, dataM, dataScore, left, right } = w;

  // The left half's last column, unwrapped top first; scoresLeftStartIdx is the row of its first cell.
  buildPeq(w, w.peq, query, queryLength);
  let status = nw(w, w.peq, W, maxNumBlocks, queryLength, target, targetLength, bestScore, false, leftHalfWidth - 1);
  if (status !== STOPPED) throw new Error('myers: no optimal path through the left half (an internal error)');
  const firstLeft = w.firstBlocks[0];
  const lastLeft = w.lastBlocks[0];
  let scoresLeftLength = (lastLeft - firstLeft + 1) * WORD_SIZE;
  for (let b = firstLeft; b <= lastLeft; b++) {
    const at = (b - firstLeft) * WORD_SIZE;
    // readBlock: dest[WORD_SIZE - 1 - i] is cell i from the bottom.
    let s = dataScore[b];
    for (let i = 0, p = dataP[b], m = dataM[b]; i < WORD_SIZE; i++, p <<= 1, m <<= 1) {
      left[at + WORD_SIZE - 1 - i] = s;
      s += (m >>> 31) - (p >>> 31);
    }
  }
  const scoresLeftStartIdx = firstLeft * WORD_SIZE;
  if (lastLeft === maxNumBlocks - 1) scoresLeftLength -= W;

  // The right half's first column from the reversed problem, unwrapped and reversed (top first in the forward rows).
  buildPeq(w, w.rPeq, rQuery, queryLength);
  status = nw(w, w.rPeq, W, maxNumBlocks, queryLength, rTarget, targetLength, bestScore, false, rightHalfWidth - 1);
  if (status !== STOPPED) throw new Error('myers: no optimal path through the right half (an internal error)');
  const firstRight = w.firstBlocks[0];
  const lastRight = w.lastBlocks[0];
  let scoresRightLength = (lastRight - firstRight + 1) * WORD_SIZE;
  for (let b = firstRight; b <= lastRight; b++) {
    const at = (lastRight - b) * WORD_SIZE;
    // readBlockReverse: dest[i] is cell i from the bottom.
    let s = dataScore[b];
    for (let i = 0, p = dataP[b], m = dataM[b]; i < WORD_SIZE; i++, p <<= 1, m <<= 1) {
      right[at + i] = s;
      s += (m >>> 31) - (p >>> 31);
    }
  }
  let scoresRightStartIdx = queryLength - (lastRight + 1) * WORD_SIZE;
  // Padding at the start of the reversed scores: skip it.
  let rightOffset = 0;
  if (scoresRightStartIdx < 0) {
    rightOffset = W;
    scoresRightStartIdx += W;
    scoresRightLength -= W;
  }

  // The row of the left column's cell which, with its lower right neighbor in the right column, sums to the best
  // score; the boundary cells (row -1, and the right column's boundary below the last row) checked last.
  const queryIdxLeftStart = Math.max(scoresLeftStartIdx, scoresRightStartIdx - 1);
  const queryIdxLeftEnd = Math.min(
    scoresLeftStartIdx + scoresLeftLength - 1,
    scoresRightStartIdx + scoresRightLength - 2,
  );
  let leftScore = -1;
  let rightScore = -1;
  let queryIdxLeftAlignment = -2; // -2: not found
  for (let queryIdx = queryIdxLeftStart; queryIdx <= queryIdxLeftEnd; queryIdx++) {
    leftScore = left[queryIdx - scoresLeftStartIdx];
    rightScore = right[rightOffset + queryIdx + 1 - scoresRightStartIdx];
    if (leftScore + rightScore === bestScore) {
      queryIdxLeftAlignment = queryIdx;
      break;
    }
  }
  if (queryIdxLeftAlignment === -2 && scoresLeftStartIdx === 0 && scoresRightStartIdx === 0) {
    leftScore = leftHalfWidth;
    rightScore = right[rightOffset];
    if (leftScore + rightScore === bestScore) queryIdxLeftAlignment = -1;
  }
  if (
    queryIdxLeftAlignment === -2 &&
    scoresLeftStartIdx + scoresLeftLength === queryLength &&
    scoresRightStartIdx + scoresRightLength === queryLength
  ) {
    leftScore = left[scoresLeftLength - 1];
    rightScore = rightHalfWidth;
    if (leftScore + rightScore === bestScore) queryIdxLeftAlignment = queryLength - 1;
  }
  if (queryIdxLeftAlignment === -2) throw new Error('myers: no optimal path across the halves (an internal error)');

  // The upper left and the lower right quarter, in order: their paths concatenate.
  const ulHeight = queryIdxLeftAlignment + 1;
  const lrHeight = queryLength - ulHeight;
  obtainAlignment(
    w,
    query,
    rQuery.subarray(lrHeight),
    ulHeight,
    target,
    rTarget.subarray(rightHalfWidth),
    leftHalfWidth,
    leftScore,
  );
  obtainAlignment(
    w,
    query.subarray(ulHeight),
    rQuery,
    lrHeight,
    target.subarray(leftHalfWidth),
    rTarget,
    rightHalfWidth,
    rightScore,
  );
}
