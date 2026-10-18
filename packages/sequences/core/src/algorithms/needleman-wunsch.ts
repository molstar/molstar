/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `needleman_wunsch` aligner: global alignment with affine gaps (Gotoh), BLOSUM62 for amino acids, a match and
// mismatch score otherwise. A gap of k letters costs `gap_open + (k - 1) * gap_extend`. Moved from @molstar/seq-core/model's
// `globalAlign` (M2b), its defaults unchanged.

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

const NEG = -1e9;
/** The first of three best: ties prefer a pair, then a gap in b. */
const argmax = (a: number, b: number, c: number) => (a >= b && a >= c ? 0 : b >= c ? 1 : 2);

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
    return globalAlign(A, B, {
      open: o.gap_open as number,
      extend: o.gap_extend as number,
      free: o.end_gaps === 'free',
      match: o.match as number,
      mismatch: o.mismatch as number,
      useMatrix,
    });
  },
};

interface Scoring {
  open: number;
  extend: number;
  free: boolean;
  match: number;
  mismatch: number;
  useMatrix: boolean;
}

/** Ties prefer a pair, then a gap in b. With free end gaps, gaps in the first and last row and column cost nothing. */
function globalAlign(a: string, b: string, s: Scoring): string {
  const n = a.length;
  const m = b.length;
  checkCells('needleman_wunsch', n, m);
  const W = m + 1;
  const free = s.free;
  // Traceback per cell: bits 0-1 M's source, 2-3 I's (consuming a), 4-5 D's (consuming b); states 0 M, 1 I, 2 D.
  const tb = new Uint8Array((n + 1) * W);
  let pM = new Float64Array(W).fill(NEG);
  let pI = new Float64Array(W).fill(NEG);
  let pD = new Float64Array(W);
  let cM = new Float64Array(W);
  let cI = new Float64Array(W);
  let cD = new Float64Array(W);
  pM[0] = 0;
  for (let j = 1; j <= m; j++) {
    pD[j] = free ? 0 : -s.open - (j - 1) * s.extend;
    tb[j] = 2 << 4;
  }
  pD[0] = NEG;
  for (let i = 1; i <= n; i++) {
    const x = a.charCodeAt(i - 1);
    // Gaps in b along the last row (trailing D) are free with free end gaps.
    const oD = free && i === n ? 0 : s.open;
    const eD = free && i === n ? 0 : s.extend;
    cM[0] = NEG;
    cD[0] = NEG;
    cI[0] = free ? 0 : -s.open - (i - 1) * s.extend;
    tb[i * W] = (i > 1 ? 1 : 0) << 2;
    for (let j = 1; j <= m; j++) {
      const y = b.charCodeAt(j - 1);
      const sc = s.useMatrix ? blosum62(x, y) : x === y ? s.match : s.mismatch;
      // Gaps in a along the last column (trailing I) are free with free end gaps.
      const oI = free && j === m ? 0 : s.open;
      const eI = free && j === m ? 0 : s.extend;
      const t = argmax(pM[j - 1], pI[j - 1], pD[j - 1]);
      cM[j] = Math.max(pM[j - 1], pI[j - 1], pD[j - 1]) + sc;
      const ti = argmax(pM[j] - oI, pI[j] - eI, pD[j] - oI);
      cI[j] = Math.max(pM[j] - oI, pI[j] - eI, pD[j] - oI);
      const td = argmax(cM[j - 1] - oD, cI[j - 1] - oD, cD[j - 1] - eD);
      cD[j] = Math.max(cM[j - 1] - oD, cI[j - 1] - oD, cD[j - 1] - eD);
      tb[i * W + j] = t | (ti << 2) | (td << 4);
    }
    // Swap the row buffers with temporaries: no array per row.
    let t = pM;
    pM = cM;
    cM = t;
    t = pI;
    pI = cI;
    cI = t;
    t = pD;
    pD = cD;
    cD = t;
  }
  let state = argmax(pM[m], pI[m], pD[m]);
  const ops: string[] = [];
  let i = n;
  let j = m;
  while (i > 0 || j > 0) {
    const cell = tb[i * W + j];
    const op = i === 0 ? 2 : j === 0 ? 1 : state;
    ops.push('MID'[op]);
    const next = (cell >> (2 * op)) & 3;
    if (op !== 2) i--;
    if (op !== 1) j--;
    state = next;
  }
  return cigarOf(ops.reverse().join(''));
}
