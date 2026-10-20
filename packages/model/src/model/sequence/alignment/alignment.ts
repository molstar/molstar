/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ryan DiRisio <rjdiris@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AlignmentOp, globalAlign } from '@molstar/core/data/util/needleman-wunsch';
import { isDebugMode } from '@molstar/core/util/debug';
import { type SubstitutionMatrix, SubstitutionMatrices } from './substitution-matrix.js';

const DefaultAlignmentOptions = {
  gapPenalty: -11,
  gapExtensionPenalty: -1,
  substMatrix: 'default' as SubstitutionMatrix | 'default',
};
export type AlignmentOptions = typeof DefaultAlignmentOptions;

/**
 * Global alignment of two sequences of residue names: `aliA` and `aliB` are the gapped rows (`-` a gap), `score` the
 * alignment's score. Without a substitution matrix equal names score 5 and different ones -3; a pair a matrix lacks
 * scores 0. A gap of k residues scores `gapPenalty + (k - 1) * gapExtensionPenalty`, except before the first pair, where
 * it scores `gapPenalty` whatever its length.
 */
export function align(seqA: ArrayLike<string>, seqB: ArrayLike<string>, options: Partial<AlignmentOptions> = {}) {
  const o = { ...DefaultAlignmentOptions, ...options };
  const substMatrix = o.substMatrix === 'default' ? undefined : SubstitutionMatrices[o.substMatrix];

  // Residue names as small integer codes, shared by both sequences, so equal names have equal codes.
  const codes = new Map<string, number>();
  const encode = (seq: ArrayLike<string>) => {
    const out = new Int32Array(seq.length);
    for (let i = 0, il = seq.length; i < il; ++i) {
      let c = codes.get(seq[i]);
      if (c === undefined) codes.set(seq[i], (c = codes.size));
      out[i] = c;
    }
    return out;
  };
  const a = encode(seqA);
  const b = encode(seqB);

  let matrix: { size: number; scores: Float64Array } | undefined;
  if (substMatrix) {
    const names = [...codes.keys()];
    const size = names.length;
    const scores = new Float64Array(size * size);
    // the pairs the DP compares, for the debug warning
    const inA = new Uint8Array(size);
    const inB = new Uint8Array(size);
    for (let i = 0; i < a.length; ++i) inA[a[i]] = 1;
    for (let j = 0; j < b.length; ++j) inB[b[j]] = 1;
    for (let p = 0; p < size; ++p) {
      for (let q = 0; q < size; ++q) {
        const val = substMatrix[names[p]]?.[names[q]];
        if (isDebugMode && val === undefined && inA[p] && inB[q]) {
          console.warn(
            `Substitution matrix has no entry for residue pair ('${names[p]}', '${names[q]}'). Non-standard residues are not supported with BLOSUM matrices.`,
          );
        }
        scores[p * size + q] = val ?? 0;
      }
    }
    matrix = { size, scores };
  }

  // The previous own DP (a number[][] S, V and H; perf-audit 04/P1) let a gap state reopen from itself and preferred
  // extending a gap on ties. When an extension costs more than an opening that is the same as extending at the opening's
  // cost with opening preferred on ties, which keeps its results identical for every option.
  const open = -o.gapPenalty;
  const extend = -o.gapExtensionPenalty;
  const { ops, score } = globalAlign(a, b, {
    matrix,
    match: 5,
    mismatch: -3,
    gapOpen: open,
    gapExtend: Math.min(open, extend),
    leadingGaps: 'flat',
    trailingGaps: 'penalized',
    gapTies: extend > open ? 'open' : 'extend',
    // as the previous DP's impossible states
    impossible: -Infinity,
  });

  const aliA: string[] = new Array(ops.length);
  const aliB: string[] = new Array(ops.length);
  for (let k = 0, i = 0, j = 0, kl = ops.length; k < kl; ++k) {
    const op = ops[k];
    aliA[k] = op === AlignmentOp.OnlyB ? '-' : seqA[i++];
    aliB[k] = op === AlignmentOp.OnlyA ? '-' : seqB[j++];
  }
  // The previous DP scored an empty `seqB` against a non-empty `seqA` as max(gapPenalty, 0) (its V matrix's column 0 was
  // left at 0, -0 kept); kept as it was.
  const emptyB = seqB.length === 0 && seqA.length > 0;
  return { aliA, aliB, score: emptyB ? (o.gapPenalty >= 0 ? o.gapPenalty : 0) : score };
}
