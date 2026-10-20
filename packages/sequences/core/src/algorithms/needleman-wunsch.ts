/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `needleman_wunsch` aligner: global alignment with affine gaps (Gotoh), BLOSUM62 for amino acids, a match and
// mismatch score otherwise. A gap of k letters costs `gap_open + (k - 1) * gap_extend`. Moved from @molstar/seq-core/model's
// `globalAlign` (M2b), its defaults unchanged. The DP is Mol*'s shared one (`@molstar/core/algorithms/needleman-wunsch`,
// also behind `@molstar/model`'s `align`): two rows of scores and one byte of traceback a cell.

import { globalAlign } from '@molstar/core/algorithms/needleman-wunsch';
import { type Aligner, type AlignerParams, alignable, cigarOf, pairType, readOptions } from './aligner.js';
import { checkCells } from './limits.js';

const ORDER = 'ARNDCQEGHILKMFPSTWYVBZX*';
// biome-ignore format: BLOSUM62's lower triangle, one row per letter in ORDER
const LOWER = [
  "4", "-1 5", "-2 0 6", "-2 -2 1 6", "0 -3 -3 -3 9", "-1 1 0 0 -3 5", "-1 0 0 2 -4 2 5", "0 -2 0 -1 -3 -2 -2 6",
  "-2 0 1 -1 -3 0 0 -2 8", "-1 -3 -3 -3 -1 -3 -3 -4 -3 4", "-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4",
  "-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5", "-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5", "-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6",
  "-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7", "1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4",
  "0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5", "-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11",
  "-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7", "0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4",
  "-2 -1 3 4 -3 0 1 -1 0 -3 -4 0 -3 -3 -2 0 -1 -4 -3 -3 4", "-1 0 0 1 -3 3 4 -2 0 -3 -3 1 -1 -3 -1 0 -1 -3 -2 -2 1 4",
  "0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 0 0 -2 -1 -1 -1 -1 -1",
  "-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 1",
];
const K = ORDER.length;
/** Code → row in ORDER; letters BLOSUM62 lacks (J, U, O) score as X. */
const INDEX = new Int8Array(128).fill(ORDER.indexOf('X'));
for (let i = 0; i < K; i++) INDEX[ORDER.charCodeAt(i)] = i;
const BLOSUM62 = new Int8Array(K * K);
for (let i = 0; i < K; i++) {
  const row = LOWER[i].split(' ');
  for (let j = 0; j <= i; j++) BLOSUM62[i * K + j] = BLOSUM62[j * K + i] = Number(row[j]);
}

/** BLOSUM62 for two upper-case letters (ASCII codes). */
export function blosum62(a: number, b: number): number {
  return BLOSUM62[INDEX[a & 127] * K + INDEX[b & 127]];
}

const params: AlignerParams = {
  matrix: {
    type: 'choice',
    choices: ['BLOSUM62'],
    nullable: true,
    default: null,
    description: 'The substitution matrix; null: BLOSUM62 for amino acids, `match`/`mismatch` for other letters.',
  },
  gap_open: { type: 'float', nullable: false, min: 0, default: 10, description: "The cost of a gap's first letter." },
  gap_extend: { type: 'float', nullable: false, min: 0, default: 1, description: 'The cost of each further letter.' },
  end_gaps: {
    type: 'choice',
    choices: ['penalized', 'free'],
    nullable: false,
    default: 'penalized',
    description: 'free: gaps before and after either sequence cost nothing (a fragment placed in a longer sequence).',
  },
  match: { type: 'float', nullable: false, default: 5, description: 'Without a matrix: the score of equal letters.' },
  mismatch: {
    type: 'float',
    nullable: false,
    default: -4,
    description: 'Without a matrix: the score of different letters.',
  },
  type: {
    type: 'choice',
    choices: ['auto', 'amino', 'nucleic'],
    nullable: false,
    default: 'auto',
    description: 'What the sequences are, for the default matrix; auto: detected (amino when both are amino acids).',
  },
};

export const needlemanWunsch: Aligner = {
  name: 'needleman_wunsch',
  description:
    'Global alignment with affine gaps (Gotoh). Memory: one byte per cell of the full matrix, refused above 2^26 cells.',
  params,
  async align(a, b, options) {
    const o = readOptions('needleman_wunsch', params, options);
    if (!alignable(a, b)) return null;
    const A = a.toUpperCase();
    const B = b.toUpperCase();
    const useMatrix = o.matrix === 'BLOSUM62' || (o.matrix === null && pairType(A, B, o.type) === 'amino');
    return align(A, B, {
      open: o.gap_open as number,
      extend: o.gap_extend as number,
      free: o.end_gaps === 'free',
      match: o.match as number,
      mismatch: o.mismatch as number,
      useMatrix,
    });
  },
};

/** Ties prefer a pair, then a gap in b. With free end gaps, gaps in the first and last row and column cost nothing. */
function align(a: string, b: string, s: Scoring): string {
  checkCells('needleman_wunsch', a.length, b.length);
  // BLOSUM62's row of each letter, or the letter's code for match and mismatch.
  const code = (seq: string) => {
    const out = new Int32Array(seq.length);
    for (let i = 0; i < seq.length; i++) out[i] = s.useMatrix ? INDEX[seq.charCodeAt(i) & 127] : seq.charCodeAt(i);
    return out;
  };
  const { ops } = globalAlign(code(a), code(b), {
    matrix: s.useMatrix ? { size: K, scores: BLOSUM62 } : undefined,
    match: s.match,
    mismatch: s.mismatch,
    gapOpen: s.open,
    gapExtend: s.extend,
    leadingGaps: s.free ? 'free' : 'penalized',
    trailingGaps: s.free ? 'free' : 'penalized',
    gapTies: 'open',
    // its own DP's sentinel, kept so that huge costs align as before
    impossible: -1e9,
  });
  let columns = '';
  for (let k = 0; k < ops.length; k++) columns += 'MID'[ops[k]];
  return cigarOf(columns);
}

interface Scoring {
  open: number;
  extend: number;
  free: boolean;
  match: number;
  mismatch: number;
  useMatrix: boolean;
}
