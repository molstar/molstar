/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Codon maps (_spec/README.md "Resolve", Maps): a protein on its gene, from exons and a strand, as the runs of its
// coding sequence: fine place `c` (residue `c / 3`) to a nucleotide. On the minus strand the coding sequence starts at
// the highest exon's end and reads downward (EC-15), so an exon junction may split a codon. And the status of a
// residue against the letters of both, through an NCBI genetic code (EC-16).

import { ImplicitLinear, type Linear, type Run, RunsLinear } from './linear.js';

/**
 * The exons as genomic pieces: an exon that ends before it starts crosses the origin of a cyclic gene of `size`
 * positions and is two pieces, `[start, size)` then `[0, end)` (EC-11); without a size it is left out.
 */
function pieces(starts: ArrayLike<number>, ends: ArrayLike<number>, size?: number): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < Math.min(starts.length, ends.length); i++) {
    if (ends[i] >= starts[i]) out.push([starts[i], ends[i]]);
    else if (size !== undefined) out.push([starts[i], size], [0, ends[i]]);
  }
  return out;
}

/**
 * The runs of a coding sequence: exons in genomic order, `strand` the way the protein runs; on the minus strand from
 * the last piece downward. `size` is a cyclic gene's length, for an exon over its origin.
 */
export function codonRuns(starts: ArrayLike<number>, ends: ArrayLike<number>, strand: '+' | '-', size?: number): Run[] {
  const list = pieces(starts, ends, size);
  if (strand === '-') list.reverse();
  const runs: Run[] = [];
  let c = 0;
  for (const [s, e] of list) {
    const length = e - s;
    if (length <= 0) continue;
    runs.push(strand === '+' ? { from: c, to: s, length, dir: 1 } : { from: c, to: e - 1, length, dir: -1 });
    c += length;
  }
  return runs;
}

/** Whether an exon crosses the origin and no size says where it splits. */
export function hasUnsplitExon(starts: ArrayLike<number>, ends: ArrayLike<number>, size?: number): boolean {
  if (size !== undefined) return false;
  for (let i = 0; i < Math.min(starts.length, ends.length); i++) if (ends[i] < starts[i]) return true;
  return false;
}

/** The cheapest encoding of a codon map's linear part: arithmetic within one piece, else its runs. */
export function codonLinear(
  starts: ArrayLike<number>,
  ends: ArrayLike<number>,
  strand: '+' | '-',
  size?: number,
): Linear {
  const runs = codonRuns(starts, ends, strand, size);
  if (runs.length === 1) {
    const r = runs[0];
    return new ImplicitLinear(r.to, r.dir, 0, r.length);
  }
  return new RunsLinear(runs);
}

/** How many nucleotides the coding sequence has. */
export function codingLength(runs: readonly Run[]): number {
  let n = 0;
  for (const r of runs) n += r.length;
  return n;
}

// ---- translation ----

const Bases = 'TCAG';

/** A genetic code: amino acids of the 64 codons in TCAG order, and the codons it starts at (read as M when first). */
export interface GeneticCode {
  readonly id: number;
  readonly aminoAcids: string;
  readonly starts: string;
}

/** The NCBI genetic codes built so far, by id: 1 (standard), 2 (vertebrate mitochondrial), 11 (bacterial, plastid). */
export const GeneticCodes: Readonly<Record<number, GeneticCode>> = {
  1: {
    id: 1,
    aminoAcids: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    starts: '---M------**--*----M---------------M----------------------------',
  },
  2: {
    id: 2,
    aminoAcids: 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
    starts: '--------------------------------MMMM---------------M------------',
  },
  11: {
    id: 11,
    aminoAcids: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
    starts: '---M------**--*----M------------MMMM---------------M------------',
  },
};

/** The code a codon map without a `table` reads with: the standard one. */
export const DefaultGeneticCode = 1;

const Standard = GeneticCodes[DefaultGeneticCode];

/** IUPAC nucleotide codes as the bases they stand for (U as T). */
const Iupac: Readonly<Record<string, string>> = {
  A: 'A',
  C: 'C',
  G: 'G',
  T: 'T',
  U: 'T',
  R: 'AG',
  Y: 'CT',
  K: 'GT',
  M: 'AC',
  S: 'CG',
  W: 'AT',
  B: 'CGT',
  D: 'AGT',
  H: 'ACT',
  V: 'ACG',
  N: 'ACGT',
};

/** The concrete codons a codon of IUPAC letters stands for (an unknown letter is N). */
function concrete(codon: string): string[] {
  let out = [''];
  for (const b of codon.toUpperCase()) {
    const set = Iupac[b] ?? Iupac.N;
    out = out.flatMap((p) => [...set].map((x) => p + x));
  }
  return out;
}

function codonIndex(c: string): number {
  return Bases.indexOf(c[0]) * 16 + Bases.indexOf(c[1]) * 4 + Bases.indexOf(c[2]);
}

/**
 * The residues a codon pairs with: the amino acids of every concrete codon its IUPAC letters stand for; a stop pairs
 * with nothing, but TGA with Sec (U) and TAG with Pyl (O); the first codon pairs with M when it is a start of the code.
 */
export function residuesOf(codon: string, first: boolean, code: GeneticCode = Standard): Set<string> {
  const out = new Set<string>();
  if (codon.length !== 3) return out;
  for (const c of concrete(codon)) {
    const i = codonIndex(c);
    const aa = code.aminoAcids[i];
    if (aa === '*') {
      if (c === 'TGA') out.add('U');
      if (c === 'TAG') out.add('O');
    } else out.add(aa);
    if (first && code.starts[i] === 'M') out.add('M');
  }
  return out;
}

/** A codon's amino acid in a code: `*` a stop, X when its letters stand for several. */
export function translateCodon(codon: string, code: GeneticCode = Standard): string {
  const all = new Set(concrete(codon).map((c) => code.aminoAcids[codonIndex(c)]));
  return all.size === 1 ? [...all][0] : 'X';
}

/** Whether a residue pairs with its codon: it is one the codon may encode; an unknown residue (X) pairs with any. */
export function codonPairs(residue: string, codon: string, first: boolean, code: GeneticCode = Standard): boolean {
  const r = residue.toUpperCase();
  if (r === 'X') return codon.length === 3;
  return residuesOf(codon, first, code).has(r);
}

const Complement: Readonly<Record<string, string>> = {
  A: 'T',
  C: 'G',
  G: 'C',
  T: 'A',
  U: 'A',
  R: 'Y',
  Y: 'R',
  K: 'M',
  M: 'K',
  S: 'S',
  W: 'W',
  B: 'V',
  V: 'B',
  D: 'H',
  H: 'D',
  N: 'N',
};

/**
 * The status of each residue of a protein against its gene's letters through a codon map's linear part: `unpaired`
 * where its triplet is not whole, `conflict` where the codon cannot encode it, else `paired`. A minus-strand codon is
 * read on the complement.
 */
export function codonStatus(
  linear: Linear,
  protein: string,
  gene: string,
  code: GeneticCode = Standard,
): ('paired' | 'conflict' | 'unpaired')[] {
  const out: ('paired' | 'conflict' | 'unpaired')[] = [];
  for (let k = 0; k < protein.length; k++) {
    let codon = '';
    for (const r of linear.forward(3 * k, 3 * k + 3))
      for (let i = 0; i < r.length; i++) {
        const b = gene[r.to + r.dir * i]?.toUpperCase() ?? 'N';
        codon += r.dir === 1 ? b : (Complement[b] ?? 'N');
      }
    if (codon.length !== 3) out.push('unpaired');
    else out.push(codonPairs(protein[k], codon, k === 0, code) ? 'paired' : 'conflict');
  }
  return out;
}
