/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Shared fixture generators of the packages' tests (RB0): seeded, by value, no network and no format: they return the
// data, and a test wraps it in its package's own types (the builder's, the model's). The genome generator already
// exists (lib/seqstar/core/_test/format/genome-fixture.ts) and is ported with the resolver.

/** mulberry32: a small seeded PRNG, uniform in [0, 1). */
export function prng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const AMINO = "ACDEFGHIKLMNPQRSTVWY";
const dna = (rand: () => number, n: number) => Array.from({ length: n }, () => "ACGT"[(rand() * 4) | 0]).join("");

export interface AlignmentFixture {
  readonly length: number;
  /** Row labels and gapped rows, all `length` long. */
  readonly labels: readonly string[];
  readonly rows: readonly string[];
  /** The rows tied to a sequence of their own (a unit), the rest being text only: `sequence` is the row without gaps. */
  readonly bound: readonly { readonly row: number; readonly sequence: string }[];
}

/** A 20,000-row protein alignment by default, two rows bound (the first and the middle one), the others text only. */
export function alignment(o: { rows?: number; length?: number; seed?: number } = {}): AlignmentFixture {
  const { rows: count = 20_000, length = 240, seed = 1 } = o;
  const rand = prng(seed);
  const root = Array.from({ length }, () => AMINO[(rand() * 20) | 0]);
  const rows = Array.from({ length: count }, () =>
    root.map((c) => (rand() < 0.08 ? "-" : rand() < 0.15 ? AMINO[(rand() * 20) | 0] : c)).join(""),
  );
  const bound = [0, count >> 1].map((row) => ({ row, sequence: rows[row].replaceAll("-", "") }));
  return { length, labels: rows.map((_, i) => `seq${i}`), rows, bound };
}

export interface GeneFixture {
  /** The genomic DNA of the plus strand. */
  readonly dna: string;
  /** The exons, `[start, stop)` on `dna`, in order. The first starts the CDS and the last holds the stop codon. */
  readonly exons: readonly (readonly [number, number])[];
  /** The coding sequence (exons joined, ATG to the stop codon) and the protein it translates to (no stop). */
  readonly cds: string;
  readonly protein: string;
}

// One codon per amino acid of the standard code, so a protein and the coding sequence of it agree.
const CODON: Record<string, string> = {
  A: "GCT", C: "TGT", D: "GAT", E: "GAA", F: "TTT", G: "GGT", H: "CAT", I: "ATT", K: "AAA", L: "CTG",
  M: "ATG", N: "AAT", P: "CCT", Q: "CAA", R: "CGT", S: "TCT", T: "ACT", V: "GTT", W: "TGG", Y: "TAT",
}; // biome-ignore format: a table

/** A DMD-sized gene by default: 79 exons over 2.2 million nucleotides coding a 3,685-residue protein. */
export function gene(o: { exons?: number; span?: number; residues?: number; seed?: number } = {}): GeneFixture {
  const { exons: count = 79, span = 2_220_000, residues = 3685, seed = 2 } = o;
  const rand = prng(seed);
  const protein = `M${Array.from({ length: residues - 1 }, () => AMINO[(rand() * 20) | 0]).join("")}`;
  const cds = `${[...protein].map((c) => CODON[c]).join("")}TAA`;
  // The CDS cut into `count` pieces of random size, laid along the span with an intron before each exon.
  const cuts = Array.from({ length: count - 1 }, () => (rand() * cds.length) | 0).sort((a, b) => a - b);
  const pieces = [0, ...cuts, cds.length].slice(0, -1).map((s, i, a) => cds.slice(s, a[i + 1] ?? cds.length));
  const intron = Math.floor((span - cds.length) / count);
  let at = 0;
  const parts: string[] = [];
  const exons = pieces.map((p): [number, number] => {
    const i = `GT${dna(rand, intron - 4)}AG`;
    parts.push(i, p);
    at += i.length;
    const e: [number, number] = [at, at + p.length];
    at += p.length;
    return e;
  });
  return { dna: parts.join(""), exons, cds, protein };
}

export interface TableFixture {
  readonly rows: number;
  readonly categories: readonly string[];
  /** Columns as typed arrays, `rows` long: a start, a stop, a value and a category index. */
  readonly start: Int32Array;
  readonly stop: Int32Array;
  readonly value: Float32Array;
  readonly category: Uint8Array;
}

/** A 3,000,000-row table by default (a variant or feature track over a genome), built once and read as typed arrays. */
export function table(o: { rows?: number; seed?: number } = {}): TableFixture {
  const { rows = 3_000_000, seed = 3 } = o;
  const rand = prng(seed);
  const start = new Int32Array(rows);
  const stop = new Int32Array(rows);
  const value = new Float32Array(rows);
  const category = new Uint8Array(rows);
  let at = 0;
  for (let i = 0; i < rows; i++) {
    at += 1 + ((rand() * 400) | 0);
    start[i] = at;
    stop[i] = at + 1 + ((rand() * 30) | 0);
    value[i] = rand();
    category[i] = (rand() * 5) | 0;
  }
  return { rows, categories: ["A", "B", "C", "D", "E"], start, stop, value, category };
}
