/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { readFileSync } from 'node:fs';
import { hashString } from '@molstar/core/data/util/hash-functions';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { ALIGNOPT, pairwiseAlignment } from '../../src/algorithms/biomsa/align/align.js';
import { extractLongestPath, noalignPair } from '../../src/algorithms/biomsa/align/noalign.js';
import { getAlignmentParameters, SEQUENCE_TYPE } from '../../src/algorithms/biomsa/align/params.js';
import { profileFromMSA } from '../../src/algorithms/biomsa/sequence/profile.js';
import { distanceKimura, makeSequence } from '../../src/algorithms/biomsa/sequence/sequence.js';
import { BioMSA } from '../../src/algorithms/biomsa/index.js';
import { GAP } from '../../src/model/alignment.js';
import {
  type Aligner,
  type AlignerOptions,
  GAP_CHAR,
  AlignerOptionsError,
  AlignerSizeError,
  Aligners,
  biomsa,
  blosum62,
  cigarFromRows,
  MAX_CELLS,
  myers,
  myersAlign,
  needlemanWunsch,
  readOptions,
} from '../../src/algorithms/index.js';
import { itWithin } from '../model/timing.js';

// 1BRS from the viewer spike's fixture (UniProtKB CC BY 4.0, wwPDB CC0), as in @molstar/seq-core/model's fixtures: entity 1 is
// barnase's mature chain (P00648 48–157), entity 2 barstar without its initiator Met (P11540 2–90) with C40A/C82A.
const spike = JSON.parse(readFileSync(new URL('../fixtures/1brs.json', import.meta.url), 'utf8'));
const P00648 = spike.barnase.sequence as string;
const P11540 = spike.barstar.sequence as string;
const ENTITY_1 =
  'AQVINTFDGVADYLQTYHKLPDNYITKSEAQALGWVASKGNLADVAPGKSIGGDIFSNREGKLPGKSGRTWREADINYTSGFRNSDRILYSSDWLIYKTTDHYQTFTKIR';
const ENTITY_2 = 'KKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDALTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS';
const DNA1 = 'ACGTACGTTTGACCA';
const DNA2 = 'ACGTACGACCA';

/** How many letters of `a` and `b` a CIGAR consumes (I consumes a, D consumes b). */
function consumed(cigar: string): [number, number] {
  let a = 0;
  let b = 0;
  for (const [, n, op] of cigar.matchAll(/(\d+)([MIDX=])/g)) {
    if (op !== 'D') a += Number(n);
    if (op !== 'I') b += Number(n);
  }
  expect(cigar).toMatch(/^(\d+[MID])+$/);
  return [a, b];
}

/** A deterministic pseudo-random sequence (a linear congruential generator) and a mutated copy of it. */
function randomPair(length: number, letters: string, seed: number): [string, string] {
  let s = seed;
  const next = () => {
    s = (Math.imul(s, 1103515245) + 12345) >>> 0;
    return s >>> 8;
  };
  let a = '';
  for (let i = 0; i < length; i++) a += letters[next() % letters.length];
  let b = '';
  for (let i = 0; i < length; i++) {
    const r = next() % 100;
    if (r < 3) continue; // deletion
    b += r < 8 ? letters[next() % letters.length] : a[i];
    if (r === 99) b += 'W'; // insertion
  }
  return [a, b];
}

/** A deterministic generator of integers, from the high bits of a linear congruential generator (its low bits cycle). */
function generator(seed: number): () => number {
  let s = seed;
  return () => {
    s = (Math.imul(s, 1103515245) + 12345) >>> 0;
    return s >>> 16;
  };
}

function randomSeq(next: () => number, length: number, letters: string): string {
  let out = '';
  for (let i = 0; i < length; i++) out += letters[next() % letters.length];
  return out;
}

function highBitsDna(length: number, seed: number): string {
  return randomSeq(generator(seed), length, 'ACGT');
}

/** A letter other than `c`. */
const unlike = (c: string | undefined) => (c === 'A' ? 'C' : 'A');

/** `insert` with its first letter unlike `after`'s and its last unlike `before`'s: it can be placed only one way. */
function fence(insert: string, before: string, after: string): string {
  return unlike(after[0]) + insert.slice(1, -1) + unlike(before.at(-1));
}

/** The longest run of M in a CIGAR, where it starts in `a` and in `b`. */
function longestMatch(cigar: string): { len: number; aStart: number; bStart: number } {
  let a = 0;
  let b = 0;
  let best = { len: 0, aStart: 0, bStart: 0 };
  for (const [, n, op] of cigar.matchAll(/(\d+)([MID])/g)) {
    const len = Number(n);
    if (op === 'M' && len > best.len) best = { len, aStart: a, bStart: b };
    if (op !== 'D') a += len;
    if (op !== 'I') b += len;
  }
  return best;
}

type Score = (x: string, y: string) => number;

/** The score of a CIGAR under affine gaps: a gap of k letters costs open + (k - 1) × extend; with `free`, none at the ends. */
function cigarScore(
  cigar: string,
  a: string,
  b: string,
  sc: Score,
  open: number,
  extend: number,
  free: boolean,
): number {
  const ops = [...cigar.matchAll(/(\d+)([MID])/g)];
  let i = 0;
  let j = 0;
  let total = 0;
  for (const [, n, op] of ops) {
    const len = Number(n);
    if (op === 'M') {
      for (let t = 0; t < len; t++) total += sc(a[i + t], b[j + t]);
    } else if (!(free && (op === 'I' ? j === 0 || j === b.length : i === 0 || i === a.length))) {
      total -= open + (len - 1) * extend;
    }
    if (op !== 'D') i += len;
    if (op !== 'I') j += len;
  }
  return total;
}

/** A plain Gotoh DP over three matrices; with `free`, gaps before and after either sequence cost nothing. */
function gotoh(a: string, b: string, sc: Score, open: number, extend: number, free: boolean): number {
  const n = a.length;
  const m = b.length;
  const NEG = Number.NEGATIVE_INFINITY;
  const grid = () => Array.from({ length: n + 1 }, () => new Array<number>(m + 1).fill(NEG));
  const M = grid();
  const I = grid(); // a letter of a against a gap
  const D = grid(); // a letter of b against a gap
  M[0][0] = 0;
  const gap = (end: boolean, k: number) => (free && end ? 0 : k);
  for (let i = 0; i <= n; i++) {
    for (let j = 0; j <= m; j++) {
      if (i > 0 && j > 0)
        M[i][j] = Math.max(M[i - 1][j - 1], I[i - 1][j - 1], D[i - 1][j - 1]) + sc(a[i - 1], b[j - 1]);
      if (i > 0) {
        const o = gap(j === 0 || j === m, open);
        const e = gap(j === 0 || j === m, extend);
        I[i][j] = Math.max(M[i - 1][j] - o, D[i - 1][j] - o, I[i - 1][j] - e);
      }
      if (j > 0) {
        const o = gap(i === 0 || i === n, open);
        const e = gap(i === 0 || i === n, extend);
        D[i][j] = Math.max(M[i][j - 1] - o, I[i][j - 1] - o, D[i][j - 1] - e);
      }
    }
  }
  return Math.max(M[n][m], I[n][m], D[n][m]);
}

/**
 * biomsa's own model (align.ts), independently: substitutions from its centred matrix (the gap extension is in the
 * centre), a gap of any length costs `open`, plus `start` when it begins a sequence and `end` when it ends one. Returns
 * the optimum and the score of a path given as biomsa's estrings (positive: letters, negative: gaps).
 */
function biomsaModel(
  a: ArrayLike<number>,
  b: ArrayLike<number>,
  S: number[][],
  open: number,
  start: number,
  end: number,
) {
  const n = a.length;
  const m = b.length;
  const NEG = Number.NEGATIVE_INFINITY;
  const grid = () => Array.from({ length: n + 1 }, () => new Array<number>(m + 1).fill(NEG));
  const M = grid();
  const V = grid(); // a letter of a against a gap, in column j
  const H = grid(); // a letter of b against a gap, in row i
  M[0][0] = 0;
  const gap = (k: number, last: number) => open + (k === 0 ? start : 0) + (k === last ? end : 0);
  for (let i = 0; i <= n; i++) {
    for (let j = 0; j <= m; j++) {
      if (i > 0 && j > 0) M[i][j] = Math.max(M[i - 1][j - 1], V[i - 1][j - 1], H[i - 1][j - 1]) + S[a[i - 1]][b[j - 1]];
      if (i > 0) V[i][j] = Math.max(M[i - 1][j] + gap(j, m), H[i - 1][j] + gap(j, m), V[i - 1][j]);
      if (j > 0) H[i][j] = Math.max(M[i][j - 1] + gap(i, n), V[i][j - 1] + gap(i, n), H[i][j - 1]);
    }
  }
  const pathScore = (ea: number[], eb: number[]) => {
    const columns = (e: number[]) => e.flatMap((x) => new Array<boolean>(Math.abs(x)).fill(x > 0));
    const ca = columns(ea);
    const cb = columns(eb);
    let i = 0;
    let j = 0;
    let total = 0;
    let prev = '';
    for (let k = 0; k < ca.length; k++) {
      const state = ca[k] && cb[k] ? 'M' : ca[k] ? 'V' : 'H';
      if (state === 'M') {
        total += S[a[i]][b[j]];
        i++;
        j++;
      } else if (state === 'V') {
        if (prev !== 'V') total += gap(j, m);
        i++;
      } else {
        if (prev !== 'H') total += gap(i, n);
        j++;
      }
      prev = state;
    }
    return i === n && j === m && ca.length === cb.length ? total : Number.NaN;
  };
  return { best: Math.max(M[n][m], V[n][m], H[n][m]), pathScore };
}

afterEach(() => vi.restoreAllMocks());

describe.each(Aligners.map((al) => [al.name, al] as [string, Aligner]))('%s', (_, aligner) => {
  it("aligns barnase's precursor to its mature chain, both ways", async () => {
    // At unit cost, pairing the chain's Ala 1 with the signal peptide's Ala 8 costs as much as one 47-letter gap:
    // `myers` finds a fragment with mode `infix` (myers.test.ts).
    const [forward, backward] = aligner.name === 'myers' ? ['7I1M40I109M', '7D1M40D109M'] : ['47I110M', '47D110M'];
    expect(await aligner.align(P00648, ENTITY_1)).toBe(forward);
    expect(await aligner.align(ENTITY_1, P00648)).toBe(backward);
  });

  it("aligns barstar to 1BRS's C40A/C82A chain without its initiator Met", async () => {
    expect(await aligner.align(P11540, ENTITY_2)).toBe('1I89M');
  });

  it('aligns the same input identically, ambiguous letters included, never drawing random numbers', async () => {
    const random = vi.spyOn(Math, 'random');
    const protein: [string, string] = ['MKVBZXLLAGRTBXXQEZ', 'MKVDEQLLGRTNXAQEZ'];
    const dna: [string, string] = ['ACGTNNRYACGTACGTACGTWSKMBDHVACGTACGTACGTACGT', 'ACGTAAGTACGTNNNNACGTACGTACGTACGT'];
    const first = [await aligner.align(...protein), await aligner.align(...dna)];
    await aligner.align(P00648, ENTITY_1);
    expect([await aligner.align(...protein), await aligner.align(...dna)]).toEqual(first);
    expect(first.every((c) => typeof c === 'string')).toBe(true);
    expect(random).not.toHaveBeenCalled();
  });

  it('gives CIGARs that consume exactly both sequences', async () => {
    const pairs: [string, string][] = [
      [P00648, ENTITY_1],
      [ENTITY_2, P11540],
      [DNA1, DNA2],
      ['a', 'ACGT'],
      ['mkv', 'mkv'],
      randomPair(300, 'ACDEFGHIKLMNPQRSTVWY', 7),
      randomPair(1700, 'ACGT', 11), // longer than 1,600: biomsa's diagonal heuristic
      // Long against short: one sequence over 1,600 letters and the other under biomsa's k-mer size of 8.
      [randomPair(1700, 'ACDEFGHIKLMNPQRSTVWY', 5)[0], 'MKV'],
      ['MKV', randomPair(1700, 'ACDEFGHIKLMNPQRSTVWY', 5)[0]],
      [randomPair(1700, 'ACGT', 3)[0], 'ACGTA'],
    ];
    for (const [a, b] of pairs) {
      expect(consumed((await aligner.align(a, b)) as string)).toEqual([a.length, b.length]);
    }
  });

  it('refuses empty input and characters other than letters', async () => {
    expect(await aligner.align('', 'ACGT')).toBeNull();
    expect(await aligner.align('ACGT', '')).toBeNull();
    expect(await aligner.align('AC-GT', 'ACGT')).toBeNull();
    expect(await aligner.align('ACGT*', 'ACGT')).toBeNull();
  });

  it("rejects unknown and invalid options with the option's name", async () => {
    await expect(aligner.align(DNA1, DNA2, { gap: 3 })).rejects.toThrow(/unknown option "gap"/);
    for (const [name, p] of Object.entries(aligner.params)) {
      if (p.min !== undefined) {
        const below = { [name]: p.min - 1 };
        await expect(aligner.align(DNA1, DNA2, below)).rejects.toThrow(`"${name}" must be at least ${p.min}`);
      }
      if (p.type === 'choice') {
        await expect(aligner.align(DNA1, DNA2, { [name]: 'rna' })).rejects.toThrow(`"${name}" must be one of`);
      }
    }
  });

  it('describes its params as data whose defaults validate', () => {
    const defaults = readOptions(aligner.name, aligner.params);
    expect(Object.keys(defaults)).toEqual(Object.keys(aligner.params));
    expect(JSON.parse(JSON.stringify(aligner.params))).toEqual(aligner.params);
  });
});

describe('biomsa', () => {
  it('changes the result with its options', async () => {
    expect(await biomsa.align(DNA1, DNA2)).toBe('6M4I5M');
    expect(await biomsa.align(DNA1, DNA2, { gap_open: 800 })).toBe('4I11M');
    expect(await biomsa.align(DNA1, DNA2, { matrix: 'BLASTZ', method: 'complete' })).toBe('6M4I5M');
  });

  it('refuses sequences of different types and a matrix for the other type', async () => {
    expect(await biomsa.align('ACGTACGTAC', 'MKWFYHEELL')).toBeNull();
    expect(await biomsa.align('ACGTACGTAC', 'MKWFYHEELL', { type: 'amino' })).toMatch(/^(\d+[MID])+$/);
    await expect(biomsa.align(DNA1, DNA2, { matrix: 'VTML240', type: 'nucleic' })).rejects.toThrow(
      /VTML240 scores amino/,
    );
    await expect(biomsa.align(DNA1, DNA2, { matrix: 'VTML240', type: 'nucleic' })).rejects.toBeInstanceOf(
      AlignerOptionsError,
    );
  });

  it('takes the type from the matrix when the type is auto', async () => {
    // ACGT letters alone are detected as nucleotides; VTML240 scores them as amino acids all the same.
    expect(await biomsa.align(DNA1, DNA2, { matrix: 'VTML240' })).toBe(
      await biomsa.align(DNA1, DNA2, { type: 'amino' }),
    );
    expect(await biomsa.align('MKWFYHEELL', 'MKWFYHEELL', { matrix: 'BLASTZ' })).toBe('10M');
    expect(await biomsa.align('ACGTACGTAC', 'MKWFYHEELL', { matrix: 'VTML240' })).toMatch(/^(\d+[MID])+$/);
  });

  it('aligns fragments of one to three letters exactly', async () => {
    for (const method of ['auto', 'complete', 'diag'] as const) {
      expect(await biomsa.align('KV', 'MKV', { method })).toBe('1D2M');
      expect(await biomsa.align('M', 'MKV', { method })).toBe('1M2D');
      expect(await biomsa.align('MKV', 'KV', { method })).toBe('1I2M');
    }
  });

  it('pairs a substring with itself: its letters all paired, on equal letters', async () => {
    const letters = 'ACDEFGHIKLMNPQRSTVWY';
    const protein = randomSeq(generator(21), 200, letters);
    const dna = randomSeq(generator(22), 120, 'ACGT');
    for (const [seq, alphabet] of [
      [protein, letters],
      [dna, 'ACGT'],
    ] as const) {
      for (const length of [3, 5, 8, 20]) {
        // Short parts sit at an end: an interior one may give two half-price end gaps for a mismatch, as biomsa intends.
        for (const start of length < 8 ? [0, seq.length - length] : [0, 37, seq.length - length]) {
          const part = seq.slice(start, start + length);
          for (const method of ['auto', 'complete'] as const) {
            const cigar = (await biomsa.align(part, seq, { method })) as string;
            const m = /^(?:(\d+)D)?(\d+)M(?:(\d+)D)?$/.exec(cigar);
            expect(m, `${alphabet.length} ${length}@${start} ${method}: ${cigar}`).not.toBeNull();
            expect(Number(m?.[2])).toBe(length);
            const from = Number(m?.[1] ?? 0);
            expect(seq.slice(from, from + length), `${length}@${start} ${method} ${cigar}`).toBe(part);
          }
        }
      }
    }
  });

  it('is optimal under its own scoring, in every end-gap mode: an independent DP (B7)', () => {
    // Mode 3 is what the diagonal method uses between diagonals: a leading gap in b cost nothing there before B7.
    const next = generator(77);
    for (let k = 0; k < 400; k++) {
      const amino = k % 2 === 0;
      const type = amino ? SEQUENCE_TYPE.PROTEIN : SEQUENCE_TYPE.NUCLEIC;
      const letters = amino ? 'ACDEFGHIKLMNPQRSTVWY' : 'ACGT';
      const a = makeSequence(randomSeq(next, 1 + (next() % 12), letters), type);
      const b = makeSequence(randomSeq(next, 1 + (next() % 12), letters), type);
      const params = getAlignmentParameters(type, {});
      const opt = [0, ALIGNOPT.DISABLE_FAVOR_START_GAP, ALIGNOPT.DISABLE_FAVOR_END_GAP, 3][(k >> 1) % 4];
      const open = params.gapOP;
      const start = opt & ALIGNOPT.DISABLE_FAVOR_START_GAP ? 0 : -open / 2;
      const end = opt & ALIGNOPT.DISABLE_FAVOR_END_GAP ? 0 : -open / 2;
      const model = biomsaModel(a.encodedSeq, b.encodedSeq, params.scoringMatrix, open, start, end);
      const result = pairwiseAlignment(a, b, params, opt);
      const label = `${a.rawSeq} ${b.rawSeq} opt ${opt}`;
      expect(result.score, label).toBe(model.best);
      expect(model.pathScore(result.estrings[0], result.estrings[1]), label).toBe(model.best);
    }
  });

  it("chains more than 32,767 diagonals (B8's lTraceBack)", () => {
    // 33,000 blocks of 10 letters; b has a T after each. The seeds are given (block k's first 8 letters, k-mer value k,
    // on diagonal -k), so the chain is 33,000 diagonals long: an Int16Array traceback broke it past 32,767, leaving a
    // leading gap of about 3 × 10^5 letters to the full DP (an AlignerSizeError).
    const blocks = 33_000;
    const next = generator(3);
    let a = '';
    let b = '';
    for (let k = 0; k < blocks; k++) {
      const block = randomSeq(next, 10, 'ACGT');
      a += block;
      b += `${block}T`;
    }
    const seeds = (step: number) => {
      const kmer = new Uint16Array(blocks);
      const kmerPos = new Int32Array(blocks);
      const winPosEnd = new Int32Array(blocks);
      const map = new Map<number, number[]>();
      for (let k = 0; k < blocks; k++) {
        kmer[k] = k;
        kmerPos[k] = step * k;
        winPosEnd[k] = step * k + 8;
        map.set(k, [k]);
      }
      const store = { kmer, kmerPos, winPos: kmerPos, winPosEnd, count: blocks };
      return [map, store, new Uint16Array(0)] as Parameters<typeof noalignPair>[3];
    };
    const [sa, sb] = [makeSequence(a, SEQUENCE_TYPE.NUCLEIC), makeSequence(b, SEQUENCE_TYPE.NUCLEIC)];
    const params = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC, {});
    const [ea, eb] = noalignPair(sa, sb, params, seeds(10), seeds(11));
    expect(eb).toEqual([b.length]);
    const gaps = ea.filter((x) => x < 0);
    expect([gaps.length, gaps.reduce((s, x) => s - x, 0)]).toEqual([blocks, blocks]);
  });

  it('keeps the seeds after a tandem repeat (B11)', async () => {
    // The minimizer store was sized for random sequence: a repeat's minimizers filled it, and the seeds after it were
    // dropped (`15000M6000D`: the last 5 kb paired with the repeat).
    const r = highBitsDna(15_000, 41);
    const [head, tail, repeat] = [r.slice(0, 10_000), r.slice(10_000), 'AC'.repeat(3000)];
    expect(await biomsa.align(head + tail, head + repeat + tail)).toBe('10000M6000D5000M');
    expect(await biomsa.align(head + repeat + tail, head + tail)).toBe('10000M6000I5000M');
  });

  it('aligns across an insertion longer than the stretch before it (B12)', async () => {
    // The diagonals were sorted by begin + diagId (2 × the position in a − the position in b): the second half's came
    // first and the first half's could not follow them, so most of the first 5 kb was paired with the insertion.
    const r = highBitsDna(16_000, 5);
    const [h1, h2] = [r.slice(0, 5000), r.slice(11_000)];
    const x = fence(r.slice(5000, 11_000), h1, h2);
    expect(await biomsa.align(h1 + h2, h1 + x + h2)).toBe('5000M6000D5000M');
    expect(await biomsa.align(h1 + x + h2, h1 + h2)).toBe('5000M6000I5000M');
  });

  /** `count` seeded insertions of 1 to 10 kb into 20 kb, each aligned both ways. */
  async function insertions(count: number) {
    // A chance k-mer match inside the insertion split it (a jump to another diagonal cost nothing), and equal letters
    // at its ends made the diagonals on either side overlap, so the second could not follow the first.
    const next = generator(2026);
    for (let k = 0; k < count; k++) {
      const base = randomSeq(next, 20_000, 'ACGT');
      const insert = randomSeq(next, 1000 + (next() % 9001), 'ACGT');
      const at = next() % 20_001;
      const long = base.slice(0, at) + insert + base.slice(at);
      for (const [a, b, op] of [
        [base, long, 'D'],
        [long, base, 'I'],
      ] as const) {
        const cigar = (await biomsa.align(a, b)) as string;
        const m = new RegExp(`^(?:(\\d+)M)?${insert.length}${op}(?:\\d+M)?$`).exec(cigar);
        expect(m, `${k}: ${insert.length} at ${at}, ${cigar.slice(0, 80)}`).not.toBeNull();
        const left = Number(m?.[1] ?? 0); // equal letters at the ends may shift the run: the rest must be the base
        expect(long.slice(0, left) + long.slice(left + insert.length) === base, `${k} ${cigar}`).toBe(true);
      }
    }
  }

  it('aligns an insertion of 1 to 10 kb anywhere in 20 kb as one run (B12), 10 insertions', () => insertions(10));

  it('aligns an insertion of 1 to 10 kb anywhere in 20 kb as one run (B12)', { tags: ['heavy'] }, () => insertions(50));

  it('aligns an insertion whose ends repeat the letters around it as one run (B12)', async () => {
    // An insertion that starts with the letters after it and ends with those before it: the diagonals overlap.
    const q = highBitsDna(16_000, 126);
    const [h1, h2] = [q.slice(0, 5000), q.slice(11_000)];
    const insert = h2.slice(0, 4) + q.slice(5004, 10_996) + h1.slice(-4);
    expect(await biomsa.align(h1 + h2, h1 + insert + h2)).toMatch(/^\d+M6000D\d+M$/);
  });

  it('breaks a near-tie toward the diagonal closer to the path (B13)', () => {
    // Two equal diagonals extend P, 10 and 30 diagonals away. The penalty subtracted P's index in the list (25, after
    // 25 one-letter diagonals that chain to nothing) from their numbers: |10 - 25| > |30 - 25| chose the farther one.
    const list: Parameters<typeof extractLongestPath>[0] = [];
    for (let k = 0; k < 25; k++) list.push({ diagId: 4 * k - 950, begin: 2 * k, end: 2 * k + 1 });
    const p = { diagId: 0, begin: 100, end: 400 };
    const near = { diagId: 10, begin: 420, end: 470 };
    const far = { diagId: 30, begin: 440, end: 490 };
    list.push(p, near, far);
    expect(extractLongestPath(list, 1000, 1000)).toEqual([p, near]);
  });

  it('places an insertion at either end, or next to it, where it is (B14)', async () => {
    // The first diagonal was chained to an empty one at (0, 0), extended over chance equal letters, and the segment
    // between them aligned with end gaps at full price; the tails were appended unaligned. Here the insertion starts
    // with h's first three letters and h is mutated at its fourth: that gave `3M6000D4997M`.
    const r = highBitsDna(11_010, 13);
    const h = r.slice(0, 5000);
    const x = h.slice(0, 3) + r.slice(5003, 11_000);
    const hx = h.slice(0, 3) + unlike(h[3]) + h.slice(4);
    expect(await biomsa.align(h, x + hx)).toBe('6000D5000M');
    expect(await biomsa.align(x + hx, h)).toBe('6000I5000M');
    const y = r.slice(5000, 10_997) + h.slice(-3); // the same at the end
    const hy = h.slice(0, -4) + unlike(h.at(-4)) + h.slice(-3);
    expect(await biomsa.align(h, hy + y)).toBe('5000M6000D');
    expect(await biomsa.align(hy + y, h)).toBe('5000M6000I');
    // Five letters from either end.
    const [s, t] = [r.slice(11_000, 11_005), r.slice(11_005)];
    const [u, v] = [fence(r.slice(5000, 11_000), h, t), fence(r.slice(5000, 11_000), s, h)];
    expect(await biomsa.align(h + t, h + u + t)).toBe('5000M6000D5M');
    expect(await biomsa.align(h + u + t, h + t)).toBe('5000M6000I5M');
    expect(await biomsa.align(s + h, s + v + h)).toBe('5M6000D5000M');
    expect(await biomsa.align(s + v + h, s + h)).toBe('5M6000I5000M');
  });

  it('places a 2 kb fragment correctly in a 70 kb sequence', async () => {
    const genome = highBitsDna(70_000, 99);
    const cigar = (await biomsa.align(genome.slice(66_000, 68_000), genome)) as string;
    const run = longestMatch(cigar);
    expect(run.len).toBeGreaterThanOrEqual(1990);
    expect(Math.abs(run.bStart - 66_000)).toBeLessThanOrEqual(10);
    expect(consumed(cigar)).toEqual([2000, 70_000]);
    // Past 65,535 letters in both sequences, with one inserted letter.
    const a = highBitsDna(400_000, 7);
    expect(await biomsa.align(a, `${a.slice(0, 200_000)}T${a.slice(200_000)}`)).toBe('200000M1D200000M');
  });

  it('aligns a short pair by the diagonal method and a long sequence against a short one', async () => {
    expect(await biomsa.align('MKV', 'MKVLAG', { method: 'diag' })).toBe('3M3D');
    expect(await biomsa.align('K', 'MKV', { method: 'diag' })).toMatch(/^(\d+[MID])+$/);
    const long = randomPair(1700, 'ACDEFGHIKLMNPQRSTVWY', 5)[0];
    expect(await biomsa.align(long.slice(900, 905), long)).toBe('900D5M795D');
  });

  it('refuses a pair whose full matrix passes the cell limit, before allocating anything', async () => {
    // 200,000 × 200,000 would take 40 billion cells: allocating them would throw a RangeError or abort the process.
    const long = 'A'.repeat(200_000);
    const error = await biomsa.align(long, long, { method: 'complete' }).catch((e) => e);
    expect(error).toBeInstanceOf(AlignerSizeError);
    expect(error).not.toBeInstanceOf(RangeError);
    expect(error.cells).toBeGreaterThan(MAX_CELLS);
    // (8191 + 1) × (8191 + 1) = 2^26 + 1 cells is over the limit.
    await expect(biomsa.align('A'.repeat(8192), 'C'.repeat(8191), { method: 'complete' })).rejects.toBeInstanceOf(
      AlignerSizeError,
    );
  });

  it('accepts a pair just under the cell limit: 8191 × 8191, its full matrix filled', { tags: ['heavy'] }, async () => {
    await expect(biomsa.align('A'.repeat(8191), 'C'.repeat(8191), { method: 'complete' })).resolves.toMatch(/M/);
  });
});

describe('the gap symbol', () => {
  it("is the model's gap code", () => {
    expect(GAP_CHAR.charCodeAt(0)).toBe(GAP);
  });

  it('is a biomsa parameter: rows, profiles and distances read it (B19)', async () => {
    const seqs = ['ACGTACGTTTGACCA', 'ACGTACGACCA', 'ACGTTTGACC'];
    const dashed = await BioMSA.align(seqs, { type: 'nucleic' });
    const dotted = await BioMSA.align(seqs, { type: 'nucleic', gapchar: '.' });
    expect(dashed.some((r) => r.includes('-'))).toBe(true);
    expect(dotted).toEqual(dashed.map((r) => r.replaceAll('-', '.')));

    const params = getAlignmentParameters(SEQUENCE_TYPE.NUCLEIC);
    const weights = dashed.map(() => 1);
    expect(profileFromMSA(dotted, weights, params, 0, '.')).toEqual(profileFromMSA(dashed, weights, params));
    expect(distanceKimura(dotted, '.')).toEqual(distanceKimura(dashed));
  });
});

describe('needleman_wunsch', () => {
  it('aligns nucleotides with match and mismatch scores', async () => {
    expect(await needlemanWunsch.align(DNA1, DNA2)).toBe('6M4I5M');
    expect(await needlemanWunsch.align('acgtacgtttgacca', DNA2)).toBe('6M4I5M');
  });

  it('changes the result with its options', async () => {
    expect(await needlemanWunsch.align('TTTTTTGGGGG', 'TTTGGGG')).toBe('3M4I4M');
    expect(await needlemanWunsch.align('TTTTTTGGGGG', 'TTTGGGG', { end_gaps: 'free' })).toBe('3I7M1I');
    expect(await needlemanWunsch.align('GATTACAGATTACA', 'GATCACAGTTACA')).toBe('8M1I5M');
    const free = { match: 1, mismatch: -1, gap_open: 0, gap_extend: 0 };
    expect(await needlemanWunsch.align('GATTACAGATTACA', 'GATCACAGTTACA', free)).toBe('3M1D1I4M1I5M');
  });

  it('scores amino acids with BLOSUM62 unless told otherwise', async () => {
    const a = 'WWWHHKTAYIAKQRQ';
    const b = 'FFYKTAYIAKQRQ';
    const blosum = await needlemanWunsch.align(a, b);
    expect(await needlemanWunsch.align(a, b, { matrix: 'BLOSUM62' })).toBe(blosum);
    expect(await needlemanWunsch.align(a, b, { type: 'nucleic', mismatch: -10 })).not.toBe(blosum);
  });
  it('refuses the cell limit too, and has no inherited option names', async () => {
    const long = 'A'.repeat(200_000);
    await expect(needlemanWunsch.align(long, long)).rejects.toBeInstanceOf(AlignerSizeError);
    await expect(needlemanWunsch.align('A'.repeat(8192), 'C'.repeat(8191))).rejects.toBeInstanceOf(AlignerSizeError);
    for (const name of ['constructor', 'toString', 'hasOwnProperty', '__proto__']) {
      const options = JSON.parse(`{"${name}": 1}`);
      for (const al of Aligners) {
        await expect(al.align(DNA1, DNA2, options), `${al.name} ${name}`).rejects.toBeInstanceOf(AlignerOptionsError);
      }
    }
  });

  // Pinned from its own DP before it moved to @molstar/core's (perf-audit 04/P1): the same CIGARs, ties included.
  it('keeps its CIGARs: seeded pairs under each option set', async () => {
    const sets: [string, string, AlignerOptions][] = [
      ['nucleic', 'ACGT', {}],
      ['nucleic, two letters', 'AC', { type: 'nucleic' }],
      ['amino, BLOSUM62', 'ACDEFGHIKLMNPQRSTVWYBZXJUOacd', { type: 'amino' }],
      ['free end gaps', 'ACGT', { end_gaps: 'free' }],
      ['free end gaps, BLOSUM62', 'ACDEFGHIKLMNPQRSTVWY', { end_gaps: 'free', matrix: 'BLOSUM62' }],
      ['fractional gaps', 'ACGT', { gap_open: 2.5, gap_extend: 0.75, match: 1.5, mismatch: -0.5 }],
      ['linear gaps', 'AC', { gap_open: 3, gap_extend: 3 }],
      ['extension over open', 'ACG', { gap_open: 1, gap_extend: 4 }],
      ['no gap cost', 'ACGT', { gap_open: 0, gap_extend: 0, match: 1, mismatch: -1 }],
    ];
    const out: Record<string, string> = {};
    for (const [k, [name, letters, options]] of sets.entries()) {
      const next = generator(31 + k);
      const cigars: string[] = [];
      for (let p = 0; p < 60; p++) {
        const a = randomSeq(next, 1 + (next() % 24), letters);
        const b = randomSeq(next, 1 + (next() % 24), letters);
        cigars.push(`${await needlemanWunsch.align(a, b, options)}`);
      }
      out[name] = cigars.join(' ');
    }
    expect(Object.fromEntries(Object.entries(out).map(([k, v]) => [k, hashString(v)]))).toEqual({
      nucleic: -1559236157,
      'nucleic, two letters': -924660635,
      'amino, BLOSUM62': -1221216645,
      'free end gaps': 443062359,
      'free end gaps, BLOSUM62': 735110485,
      'fractional gaps': 1043541977,
      'linear gaps': -690344163,
      'extension over open': 853399335,
      'no gap cost': -1943623818,
    });
  });

  itWithin('is optimal: 500 seeded pairs score what an independent Gotoh DP finds', (within) =>
    within(2000, async () => {
      const next = generator(2024);
      for (let k = 0; k < 500; k++) {
        const amino = k % 2 === 0;
        const letters = amino ? 'ACDEFGHIKLMNPQRSTVWY' : 'ACGT';
        const a = randomSeq(next, 1 + (next() % 30), letters);
        const b = randomSeq(next, 1 + (next() % 30), letters);
        const free = k % 5 === 0;
        const options = {
          type: amino ? 'amino' : 'nucleic',
          end_gaps: free ? 'free' : 'penalized',
          gap_open: 1 + (next() % 12),
          gap_extend: next() % 3,
        };
        const cigar = (await needlemanWunsch.align(a, b, options)) as string;
        const sc = (x: string, y: string) => (amino ? blosum62(x.charCodeAt(0), y.charCodeAt(0)) : x === y ? 5 : -4);
        const expected = gotoh(a, b, sc, options.gap_open, options.gap_extend, free);
        expect(cigarScore(cigar, a, b, sc, options.gap_open, options.gap_extend, free), `${a} ${b} ${cigar}`).toBe(
          expected,
        );
      }
    }),
  );
});

describe('myersAlign', () => {
  it('validates its input', () => {
    expect(() => myersAlign('ACÉT', 'ACIT')).toThrow(RangeError);
    expect(() => myersAlign('AC-T', 'ACGT')).toThrow(RangeError);
    expect(() => myersAlign('ACGT', 'ACGT', 'global', -1)).toThrow(RangeError);
    expect(() => myersAlign('ACGT', 'ACGT', 'global', 1.5)).toThrow(RangeError);
    expect(() => myersAlign('ACGT', 'ACGT', 'global', null, ['N'])).toThrow(RangeError);
    expect(() => myersAlign('ACGT', 'ACGT', 'global', null, ['ÉI'])).toThrow(RangeError);
    expect(myersAlign('ACGT', 'ACGT')?.cigar).toBe('4M');
  });

  it('has no cell limit: linear memory', async () => {
    const seq = 'ACGT'.repeat(2300); // 9,200 × 9,200 cells, over MAX_CELLS
    expect(await myers.align(seq, seq)).toBe('9200M');
  });
});

describe('error classes', () => {
  it('separates option errors from size errors', async () => {
    const options = await needlemanWunsch.align(DNA1, DNA2, { gap_open: -1 }).catch((e) => e);
    expect(options).toBeInstanceOf(AlignerOptionsError);
    expect(options).toBeInstanceOf(RangeError);
    const size = new AlignerSizeError('x', MAX_CELLS, 1);
    expect(size).not.toBeInstanceOf(AlignerOptionsError);
    expect(size).not.toBeInstanceOf(RangeError);
  });
});

describe('cigarFromRows', () => {
  it('reads gapped rows, skipping columns gapped in both', () => {
    expect(cigarFromRows('MK-V--LL', '-KAV-DL-')).toBe('1I1M1D1M1D1M1I');
    expect(() => cigarFromRows('AB', 'A')).toThrow(RangeError);
  });
});
