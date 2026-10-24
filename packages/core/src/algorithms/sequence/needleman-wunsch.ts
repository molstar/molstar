/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Global pairwise alignment with affine gaps (Gotoh's Needleman–Wunsch) over integer-coded sequences: the one DP
// behind `@molstar/model`'s `align` (SEQRES mapping, superposition) and Seq*'s `needleman_wunsch` aligner. Scores are
// kept for two rows only; the traceback is one byte per cell of the full (n + 1) × (m + 1) matrix (2 bits for each of
// the three states' sources). A caller bounds the size (Seq*'s `MAX_CELLS`) where it needs to.

/** A column of an alignment: a pair, a letter of `a` only (a gap in `b`), a letter of `b` only (a gap in `a`). */
export const AlignmentOp = { Pair: 0, OnlyA: 1, OnlyB: 2 } as const;
export type AlignmentOp = (typeof AlignmentOp)[keyof typeof AlignmentOp];

export interface GlobalAlignScoring {
  /**
   * Substitution scores, row-major `size × size`, indexed by the codes of `a` (row) and `b` (column); every code must
   * be below `size`. Without it, equal codes score `match` and different ones `mismatch`.
   */
  readonly matrix?: { readonly size: number; readonly scores: ArrayLike<number> };
  readonly match?: number;
  readonly mismatch?: number;
  /** The cost of a gap's first letter (subtracted from the score). */
  readonly gapOpen: number;
  /** The cost of each further letter of a gap. */
  readonly gapExtend: number;
  /**
   * Gaps before the first pair: `penalized`, `gapOpen + (k - 1) * gapExtend` as inside; `flat`, `gapOpen` whatever
   * their length; `free`, nothing.
   */
  readonly leadingGaps: 'penalized' | 'flat' | 'free';
  /**
   * Gaps after the last pair: `penalized` as inside, or `free`. A gap against an empty sequence is leading: it costs
   * what `leadingGaps` says.
   */
  readonly trailingGaps: 'penalized' | 'free';
  /**
   * The source a gap state takes among equally good ones. `open`: a gap in `b` opens from a pair, then extends, then
   * opens from a gap in `a`; a gap in `a` opens from a pair, then from a gap in `b`, then extends (the old Seq* argmax,
   * which the pinned results follow). `extend`: extending, then opening from a pair, then from the other gap. A pair's
   * source is always a pair, then a gap in `b`, then a gap in `a`; so is the final state.
   */
  readonly gapTies: 'open' | 'extend';
  /**
   * The score of a state no alignment reaches (the DP's border), -Infinity by default. A finite one (Seq*'s -1e9) is
   * reachable when gap costs or scores come near it.
   */
  readonly impossible?: number;
}

export interface GlobalAlignment {
  /** The columns, first to last (`AlignmentOp`); it covers both sequences. */
  readonly ops: Uint8Array;
  /** The alignment's score: pairs' scores less gap costs. */
  readonly score: number;
}

/**
 * The best global alignment of `a` and `b` under `scoring`; deterministic, ties broken as `gapTies` says. Memory is one
 * byte per cell of (n + 1) × (m + 1) plus O(m); time O(n × m).
 */
export function globalAlign(a: ArrayLike<number>, b: ArrayLike<number>, scoring: GlobalAlignScoring): GlobalAlignment {
  const n = a.length;
  const m = b.length;
  const W = m + 1;
  const open = scoring.gapOpen;
  const extend = scoring.gapExtend;
  const lead = scoring.leadingGaps;
  const freeTrail = scoring.trailingGaps === 'free';
  const extendFirst = scoring.gapTies === 'extend';
  const matrix = scoring.matrix?.scores;
  const K = scoring.matrix?.size ?? 0;
  const match = scoring.match ?? 0;
  const mismatch = scoring.mismatch ?? 0;
  const NEG = scoring.impossible ?? -Infinity;
  const leading = (k: number) => (lead === 'free' ? 0 : lead === 'flat' ? -open : -open - (k - 1) * extend);

  // Traceback per cell: bits 0-1 M's source, 2-3 I's (consuming a), 4-5 D's (consuming b); states 0 M, 1 I, 2 D.
  // Row 0 and column 0 are never read: there the op is forced.
  const tb = new Uint8Array((n + 1) * W);
  let pM = new Float64Array(W).fill(NEG);
  let pI = new Float64Array(W).fill(NEG);
  let pD = new Float64Array(W);
  let cM = new Float64Array(W);
  let cI = new Float64Array(W);
  let cD = new Float64Array(W);
  pM[0] = 0;
  pD[0] = NEG;
  for (let j = 1; j <= m; j++) pD[j] = leading(j);

  for (let i = 1; i <= n; i++) {
    const x = a[i - 1];
    const row = x * K;
    // Gaps in b along the last row (trailing D) are free with free trailing gaps.
    const oD = freeTrail && i === n ? 0 : open;
    const eD = freeTrail && i === n ? 0 : extend;
    cM[0] = NEG;
    cD[0] = NEG;
    cI[0] = leading(i);
    const o = i * W;
    for (let j = 1; j <= m; j++) {
      const y = b[j - 1];
      const sc = matrix !== undefined ? matrix[row + y] : x === y ? match : mismatch;
      // Gaps in a along the last column (trailing I) are free with free trailing gaps.
      const oI = freeTrail && j === m ? 0 : open;
      const eI = freeTrail && j === m ? 0 : extend;

      // Each state takes the first best of its three sources in tie order. Branch-free on purpose: the comparisons
      // are data-dependent, and as branches their mispredictions cost about 40 % of the time (perf-audit 04/P1).
      // The source is 2 when the last in order wins strictly, else 1 or 0 by one comparison; `gapTies: 'extend'` only
      // makes the extension's comparison non-strict (I: I, M, D; D: D, M, I), a branch on a constant.

      // M: from (i - 1, j - 1); ties M, I, D.
      const m0 = pM[j - 1];
      const m1 = pI[j - 1];
      const m2 = pD[j - 1];
      const m01 = Math.max(m0, m1);
      cM[j] = Math.max(m01, m2) + sc;
      const mD = Number(m2 > m01);
      const t = 2 * mD + Number(m1 > m0) * (1 - mD);

      // I: from (i - 1, j).
      const im = pM[j] - oI;
      const ii = pI[j] - eI;
      const id = pD[j] - oI;
      const i01 = Math.max(im, ii);
      cI[j] = Math.max(i01, id);
      const iD = Number(id > i01);
      const ti = 2 * iD + Number(extendFirst ? ii >= im : ii > im) * (1 - iD);

      // D: from (i, j - 1).
      const dm = cM[j - 1] - oD;
      const di = cI[j - 1] - oD;
      const dd = cD[j - 1] - eD;
      const d01 = Math.max(dm, di);
      cD[j] = Math.max(d01, dd);
      const dD = Number(extendFirst ? dd >= d01 : dd > d01);
      const td = 2 * dD + Number(di > dm) * (1 - dD);

      tb[o + j] = t | (ti << 2) | (td << 4);
    }
    // Swap the row buffers: no array per row.
    let s = pM;
    pM = cM;
    cM = s;
    s = pI;
    pI = cI;
    cI = s;
    s = pD;
    pD = cD;
    cD = s;
  }

  // The final state: the best at (n, m); ties M, I, D.
  let state = 0;
  let best = pM[m];
  if (pI[m] > best) {
    best = pI[m];
    state = 1;
  }
  if (pD[m] > best) state = 2;
  // `Math.max` as the cells' scores: +0 over -0
  const score = Math.max(pM[m], pI[m], pD[m]);

  const ops = new Uint8Array(n + m);
  let k = n + m;
  let i = n;
  let j = m;
  while (i > 0 || j > 0) {
    const op = i === 0 ? 2 : j === 0 ? 1 : state;
    ops[--k] = op;
    const next = (tb[i * W + j] >> (2 * op)) & 3;
    if (op !== 2) i--;
    if (op !== 1) j--;
    state = next;
  }
  return { ops: ops.subarray(k), score };
}
