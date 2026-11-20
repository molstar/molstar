/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { DerivedData } from './derived-data.js';

// The tRNA dataset: PDB 1EHZ's chain A (yeast tRNA-Phe, 76 nucleotides, 14 modified residues written as their parent
// letters, their components the sequence's `modified` table, their codes also as ranges) and its 27 base pairs, the entry's own table; Rfam RF00005's seed (954 rows,
// 118 columns, Stockholm with `SS_cons`, as published); and what is computed from them (`derive.ts`). The molecule and the
// seed are read by location, as a published file is; the pairs of the entry are held here (a `pairs` node is inline
// rows only) and compared with the committed `1EHZ-pairs.seq` by the fixtures test.

export const MoleculeFile = 'fixtures/trna/1EHZ.seq';
export const PairsFile = 'fixtures/trna/1EHZ-pairs.seq';
export const SeedFile = 'fixtures/trna/RF00005.sto';
/** The molecule's system in its file. */
export const SystemName = '1EHZ-A';
/** 1EHZ's chain A, parent letters (compared with the file's letters by the fixtures test). */
export const TrnaLetters = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA';
export const TrnaLength = TrnaLetters.length;

/**
 * 1EHZ's modified residues (`pdbx_struct_mod_residue`, 0-based): the component at a position, its letter the parent's
 * (14 positions, 11 components; compared with the file's ranges by the fixtures test).
 */
export const ModifiedResidues: readonly (readonly [position: number, component: string])[] = [
  [9, '2MG'],
  [15, 'H2U'],
  [16, 'H2U'],
  [25, 'M2G'],
  [31, 'OMC'],
  [33, 'OMG'],
  [36, 'YYG'],
  [38, 'PSU'],
  [39, '5MC'],
  [45, '7MG'],
  [48, '5MC'],
  [53, '5MU'],
  [54, 'PSU'],
  [57, '1MA'],
];

/**
 * Their compounds: each code's name (the entry's `chem_comp.name` in lower case; YYG's its common name, wybutosine, for
 * the entry's systematic one) and parent.
 */
export const ModifiedCompounds: Readonly<Record<string, { readonly name: string; readonly parent: string }>> = {
  '2MG': { name: "2N-methylguanosine-5'-monophosphate", parent: 'G' },
  H2U: { name: "5,6-dihydrouridine-5'-monophosphate", parent: 'U' },
  M2G: { name: "N2-dimethylguanosine-5'-monophosphate", parent: 'G' },
  OMC: { name: "O2'-methylycytidine-5'-monophosphate", parent: 'C' },
  OMG: { name: "O2'-methylguanosine-5'-monophosphate", parent: 'G' },
  YYG: { name: 'wybutosine', parent: 'G' },
  PSU: { name: "pseudouridine-5'-monophosphate", parent: 'U' },
  '5MC': { name: "5-methylcytidine-5'-monophosphate", parent: 'C' },
  '7MG': { name: "7N-methyl-8-hydroguanosine-5'-monophosphate", parent: 'G' },
  '5MU': { name: "5-methyluridine 5'-monophosphate", parent: 'U' },
  '1MA': { name: "6-hydro-1-methyladenosine-5'-monophosphate", parent: 'A' },
};

export type Kind = 'canonical' | 'wobble' | 'non-canonical';
const KindOf = { C: 'canonical', W: 'wobble', N: 'non-canonical' } as const;

/** The entry's base pairs: the first base, its partner (0-based, first before partner), the kind, whether tertiary. */
export const Pairs: readonly (readonly [i: number, j: number, kind: 'C' | 'W' | 'N', tertiary: 0 | 1])[] = [
  [0, 71, 'C', 0],
  [1, 70, 'C', 0],
  [2, 69, 'C', 0],
  [3, 68, 'W', 0],
  [4, 67, 'C', 0],
  [5, 66, 'C', 0],
  [6, 65, 'C', 0],
  [7, 13, 'N', 1],
  [9, 24, 'C', 0],
  [10, 23, 'C', 0],
  [11, 22, 'C', 0],
  [12, 21, 'C', 0],
  [14, 47, 'N', 1],
  [17, 54, 'N', 1],
  [18, 55, 'C', 1],
  [25, 43, 'N', 1],
  [26, 42, 'C', 0],
  [27, 41, 'C', 0],
  [28, 40, 'C', 0],
  [29, 39, 'C', 0],
  [31, 37, 'N', 0],
  [48, 64, 'C', 0],
  [49, 63, 'C', 0],
  [50, 62, 'C', 0],
  [51, 61, 'C', 0],
  [52, 60, 'C', 0],
  [53, 57, 'N', 1],
];

/** A `pairs` node's rows over a list of pairs, 0-based, an id `i-j` 1-based. */
export function pairRows(pairs: readonly (readonly [number, number, Kind, boolean?])[]) {
  return {
    // `tertiary` rides along as a column of the rows (a filter reads it: `(item :tertiary yes)`)
    ...(pairs.some((p) => p[3] !== undefined) ? { tertiary: pairs.map((p) => (p[3] ? 'yes' : 'no')) } : {}),
    a_start: pairs.map((p) => p[0]),
    a_end: pairs.map((p) => p[0] + 1),
    b_start: pairs.map((p) => p[1]),
    b_end: pairs.map((p) => p[1] + 1),
    type: pairs.map((p) => p[2]),
    id: pairs.map((p) => `${p[0] + 1}-${p[1] + 1}`),
  };
}

/** 1EHZ's base pairs as a pairs table. */
export const ehzPairs = () => pairRows(Pairs.map(([i, j, k, t]) => [i, j, KindOf[k], t === 1]));

/** The seed row 1EHZ is placed on (computed): its `SS_cons` pairs, nested only. */
export const seedRowPairs = () => pairRows(DerivedData.rowPairs.map(([i, j]) => [i, j, 'canonical']));

/**
 * The covariation of 1EHZ's positions (computed on them from the seed's columns), symmetric: the grid's values where
 * the evidence is, NaN (under the floor, or aligned to a gap) elsewhere.
 */
export function covariationGrid(): Float32Array {
  const n = TrnaLength;
  const v = new Float32Array(n * n).fill(Number.NaN);
  for (const [i, j, bits] of DerivedData.covariation) v[i * n + j] = v[j * n + i] = bits;
  return v;
}

/** The scale variant's size: a 4,000 × 4,000 grid (performance.md, Scale). */
export const ScaleSize = 4000;

/** The letters of a generated RNA of `n` nucleotides from a seed. */
export function scaleLetters(n: number): string {
  const r = new PCG(4000);
  let s = '';
  for (let i = 0; i < n; i++) s += 'ACGU'[r.int() % 4];
  return s;
}

/**
 * The scale variant's values, `n × n` and symmetric: about one cell in eight holds a small value, a tenth of `n`
 * helices of high ones, NaN elsewhere (the density of a probability matrix's signal); from a seed.
 */
export function scaleGrid(n: number): Float32Array {
  const r = new PCG(4001);
  const v = new Float32Array(n * n).fill(Number.NaN);
  for (let i = 0; i < n; i++)
    for (let j = i + 4; j < n; j++) {
      if (r.int() % 8 !== 0) continue;
      v[i * n + j] = v[j * n + i] = (r.int() % 1000) / 4000;
    }
  // helices: runs of high values along an anti-diagonal, the structure a zoomed-out plot shows
  for (let k = 0; k < n / 10; k++) {
    const i = r.int() % n;
    const j = r.int() % n;
    const len = 6 + (r.int() % 14);
    for (let h = 0; h < len; h++) {
      const a = i + h;
      const b = j - h;
      if (a >= 0 && a < n && b >= 0 && b < n && a !== b) v[a * n + b] = v[b * n + a] = 0.95;
    }
  }
  return v;
}
