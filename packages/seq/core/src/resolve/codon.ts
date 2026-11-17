/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Codon maps (_spec/README.md "Resolve", Maps): a protein on its gene, from exons and a strand, as the runs of its
// coding sequence: fine place `c` (residue `c / 3`) to a nucleotide. On the minus strand the coding sequence starts at
// the highest exon's end and reads downward (EC-15), so an exon junction may split a codon. And the status of a
// residue against the letters of both, through the standard genetic code (EC-16).

import { ImplicitLinear, type Linear, type Run, RunsLinear } from './linear.js';

/** The runs of a coding sequence: exons in genomic order, `strand` the way the protein runs. */
export function codonRuns(starts: ArrayLike<number>, ends: ArrayLike<number>, strand: '+' | '-'): Run[] {
  const n = Math.min(starts.length, ends.length);
  const order = Array.from({ length: n }, (_, i) => i);
  if (strand === '-') order.reverse();
  const runs: Run[] = [];
  let c = 0;
  for (const i of order) {
    const length = ends[i] - starts[i];
    if (length <= 0) continue;
    runs.push(
      strand === '+' ? { from: c, to: starts[i], length, dir: 1 } : { from: c, to: ends[i] - 1, length, dir: -1 },
    );
    c += length;
  }
  return runs;
}

/** The cheapest encoding of a codon map's linear part: arithmetic within one exon, else its runs. */
export function codonLinear(starts: ArrayLike<number>, ends: ArrayLike<number>, strand: '+' | '-'): Linear {
  const runs = codonRuns(starts, ends, strand);
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
// NCBI table 1, codons in TCAG order
const Standard = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
/** Codons the standard code starts at (read as M when first). */
const Starts = new Set(['TTG', 'CTG', 'ATG']);
const Complement: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A' };

/** A codon's amino acid in the standard code: T and U alike, an ambiguous base gives X, a stop `*`. */
export function translateCodon(codon: string): string {
  const c = codon.toUpperCase().replace(/U/g, 'T');
  let i = 0;
  for (const b of c) {
    const k = Bases.indexOf(b);
    if (k < 0) return 'X';
    i = i * 4 + k;
  }
  return c.length === 3 ? Standard[i] : 'X';
}

/**
 * Whether a residue pairs with its codon: equal letters; X (an ambiguous codon) pairs with any; the first codon is
 * read as M when it is a start; TGA pairs with Sec (U), TAG with Pyl (O). A stop pairs with nothing.
 */
export function codonPairs(residue: string, codon: string, first: boolean): boolean {
  const r = residue.toUpperCase();
  const c = codon.toUpperCase().replace(/U/g, 'T');
  if (first && r === 'M' && Starts.has(c)) return true;
  const aa = translateCodon(c);
  if (aa === 'X') return true;
  if (aa === '*') return (c === 'TGA' && r === 'U') || (c === 'TAG' && r === 'O');
  return aa === r;
}

/**
 * The status of each residue of a protein against its gene's letters through a codon map's linear part: `unpaired`
 * where its triplet is not whole, `conflict` where the codon does not encode it, else `paired`. A minus-strand codon
 * is read on the complement.
 */
export function codonStatus(linear: Linear, protein: string, gene: string): ('paired' | 'conflict' | 'unpaired')[] {
  const out: ('paired' | 'conflict' | 'unpaired')[] = [];
  for (let k = 0; k < protein.length; k++) {
    let codon = '';
    for (const r of linear.forward(3 * k, 3 * k + 3))
      for (let i = 0; i < r.length; i++) {
        const b = gene[r.to + r.dir * i]?.toUpperCase() ?? 'N';
        codon += r.dir === 1 ? b : (Complement[b.replace('U', 'T')] ?? 'N');
      }
    if (codon.length !== 3) out.push('unpaired');
    else out.push(codonPairs(protein[k], codon, k === 0) ? 'paired' : 'conflict');
  }
  return out;
}
