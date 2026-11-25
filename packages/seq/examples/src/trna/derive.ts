/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AlignmentOp, globalAlign } from '@molstar/core/algorithms/sequence/needleman-wunsch';
import { readStockholm } from '@molstar/seqx-core/data/alignment-readers';

// What the tRNA dataset computes from its committed fixtures: 1EHZ's place in the family (a pairwise
// alignment of its parent letters to the seed row of highest identity: `globalAlign`), the covariation between its
// positions (mutual information over the seed's columns they align to, less the Miller–Madow correction for small
// counts), and the best row's own letters and pairs (the seed's `SS_cons` projected onto it). Generated data on
// the views' side: `derived-data.ts` holds the result, written by the fixtures test (`SEQ_RECORDS=update`) and compared
// with this code by it, so a browser never reads the 141 KB seed to build a view.

const Code: Record<string, number> = { A: 0, C: 1, G: 2, U: 3, T: 3 };
const code = (c: string) => Code[c.toUpperCase()] ?? 4;

/** A pair of 1EHZ's positions and its covariation (bits), 0-based, `i < j`. */
export type CovariationEntry = readonly [i: number, j: number, bits: number];

export interface Derived {
  /** The seed row 1EHZ is placed on: its label, its aligned letters. */
  readonly row: string;
  /** Its identity with 1EHZ over 1EHZ's length. */
  readonly identity: number;
  /** The row's own letters, no gaps (U for T). */
  readonly rowLetters: string;
  /** The `SS_cons` pairs on those letters (0-based, `i < j`). */
  readonly rowPairs: readonly (readonly [number, number])[];
  /** Each of 1EHZ's positions' column of the seed (-1: aligned to a gap). */
  readonly columns: readonly number[];
  /** The seed's `SS_cons` pairs on its columns (0-based, `i < j`): the consensus arcs over the alignment. */
  readonly seedPairs: readonly (readonly [number, number])[];
  /** The covariation at each pair of 1EHZ's positions with at least `MinRows` rows and `MinBits` bits, 3 decimals. */
  readonly covariation: readonly CovariationEntry[];
}

export const MinRows = 30;
export const MinBits = 0.05;
export const Decimals = 3;

/**
 * The mutual information (bits) of two columns of codes (A, C, G, U below 4; other codes left out) with the
 * Miller–Madow correction, `NaN` under `MinRows` rows, kept in `[0, min(H(a), H(b))]` (the correction can pass the entropy).
 */
export function mutualInformation(a: Uint8Array, b: Uint8Array): number {
  const joint = new Float64Array(16);
  let n = 0;
  for (let k = 0; k < a.length; k++) {
    if (a[k] > 3 || b[k] > 3) continue;
    joint[a[k] * 4 + b[k]]++;
    n++;
  }
  if (n < MinRows) return Number.NaN;
  const pa = new Float64Array(4);
  const pb = new Float64Array(4);
  for (let x = 0; x < 4; x++)
    for (let y = 0; y < 4; y++) {
      pa[x] += joint[x * 4 + y];
      pb[y] += joint[x * 4 + y];
    }
  let mi = 0;
  for (let x = 0; x < 4; x++)
    for (let y = 0; y < 4; y++) {
      const c = joint[x * 4 + y];
      if (c > 0) mi += (c / n) * Math.log2((c * n) / (pa[x] * pb[y]));
    }
  // Miller–Madow: the bias of the plug-in estimate, (observed cells − observed rows − observed columns + 1) / 2n ln 2
  const used = (v: Float64Array) => v.reduce((s, x) => s + (x > 0 ? 1 : 0), 0);
  const cells = joint.reduce((s, x) => s + (x > 0 ? 1 : 0), 0);
  const bias = (cells - used(pa) - used(pb) + 1) / (2 * n * Math.LN2);
  const entropy = (v: Float64Array) => v.reduce((h, c) => (c > 0 ? h - (c / n) * Math.log2(c / n) : h), 0);
  return Math.max(0, Math.min(mi - bias, entropy(pa), entropy(pb)));
}

/** The place, the best row, the covariation: from the seed's text and 1EHZ's letters. */
export function derive(seed: string, letters: string): Derived {
  const read = readStockholm(seed);
  const target = Uint8Array.from(letters, code);
  const scoring = {
    match: 2,
    mismatch: -1,
    gapOpen: 3,
    gapExtend: 1,
    leadingGaps: 'penalized',
    trailingGaps: 'penalized',
    gapTies: 'open',
  } as const;
  // the row of highest identity: matches over 1EHZ's length, the first on a tie
  let best = -1;
  let bestIdentity = -1;
  let bestOps: Uint8Array | undefined;
  read.rows.forEach((row, r) => {
    const own = Uint8Array.from(row.letters.replaceAll(/[-.]/g, ''), code);
    const ops = globalAlign(target, own, scoring).ops;
    let i = 0;
    let j = 0;
    let same = 0;
    for (const op of ops) {
      if (op === AlignmentOp.Pair) {
        if (target[i] === own[j] && target[i] < 4) same++;
        i++;
        j++;
      } else if (op === AlignmentOp.OnlyA) i++;
      else j++;
    }
    if (same / letters.length > bestIdentity) {
      bestIdentity = same / letters.length;
      best = r;
      bestOps = ops;
    }
  });
  const row = read.rows[best];
  // the row's residue index of each column, and each 1EHZ position's column through the alignment
  const columnOf: number[] = [];
  for (let c = 0; c < row.letters.length; c++) if (row.letters[c] !== '-' && row.letters[c] !== '.') columnOf.push(c);
  const columns = new Array<number>(letters.length).fill(-1);
  let i = 0;
  let j = 0;
  for (const op of bestOps!) {
    if (op === AlignmentOp.Pair) columns[i++] = columnOf[j++];
    else if (op === AlignmentOp.OnlyA) i++;
    else j++;
  }
  const rowLetters = row.letters.replaceAll(/[-.]/g, '').toUpperCase().replaceAll('T', 'U');
  const residueOf = new Map(columnOf.map((c, k) => [c, k] as const));
  const rowPairs = (read.pairs ?? [])
    .filter(([a, b]) => residueOf.has(a) && residueOf.has(b))
    .map(([a, b]) => [residueOf.get(a)!, residueOf.get(b)!] as const);
  // the covariation of 1EHZ's positions over the columns they align to
  const profile = (c: number) => Uint8Array.from(read.rows, (r) => code(r.letters[c]));
  const cache = new Map<number, Uint8Array>();
  const column = (c: number) => cache.get(c) ?? (cache.set(c, profile(c)), cache.get(c)!);
  const covariation: CovariationEntry[] = [];
  const scale = 10 ** Decimals;
  for (let a = 0; a < letters.length; a++) {
    if (columns[a] < 0) continue;
    for (let b = a + 1; b < letters.length; b++) {
      if (columns[b] < 0) continue;
      const mi = mutualInformation(column(columns[a]), column(columns[b]));
      if (mi >= MinBits) covariation.push([a, b, Math.round(mi * scale) / scale]);
    }
  }
  return {
    row: row.label,
    identity: Math.round(bestIdentity * 1000) / 1000,
    rowLetters,
    rowPairs,
    columns,
    seedPairs: (read.pairs ?? []).map(([a, b]) => [a, b] as const),
    covariation,
  };
}

/** The derived data as the module `derived-data.ts` (compact: a row an entry's worth). */
export function derivedModule(d: Derived): string {
  const arr = (v: readonly unknown[]) => `[${v.map((x) => JSON.stringify(x)).join(',')}]`;
  return `/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Derived } from './derive.js';

// Generated by \`derive.ts\` from the committed fixtures (1EHZ and the Rfam RF00005 seed), compared with it by the
// fixtures test; \`SEQ_RECORDS=update\` writes it again. Do not edit.

export const DerivedData: Derived = {
  row: ${JSON.stringify(d.row)},
  identity: ${d.identity},
  rowLetters: ${JSON.stringify(d.rowLetters)},
  rowPairs: ${arr(d.rowPairs)},
  columns: ${arr(d.columns)},
  seedPairs: ${arr(d.seedPairs)},
  covariation: ${arr(d.covariation)},
};
`;
}
