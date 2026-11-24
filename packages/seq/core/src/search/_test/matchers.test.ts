/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The built-in matchers against a naive oracle on seeded letters and patterns (both strands, a circle, overlaps, chunks),
// and each hit rule by name (_spec/README.md "Search", Built-in matchers, Hits).

import { describe, expect, it } from 'vitest';
import { DegenerateMatcher, ExactMatcher, RegexMatcher } from '../index.js';
import { oracle, randomLetters, runSearch, seeded, textTarget } from './run.js';

const Iupac: Record<string, string> = {
  A: 'A',
  C: 'C',
  G: 'G',
  T: 'T',
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
const inside = (seq: string, pat: string) => [...(Iupac[seq] ?? '?')].every((b) => (Iupac[pat] ?? '').includes(b));

const exactAt = (pattern: string) => (text: string, s: number) =>
  text.startsWith(pattern.toUpperCase().replace(/U/g, 'T'), s) ? pattern.length : 0;
const degenerateAt = (pattern: string) => (text: string, s: number) => {
  const p = pattern.toUpperCase();
  if (s + p.length > text.length) return 0;
  for (let i = 0; i < p.length; i++) if (!inside(text[s + i], p[i])) return 0;
  return p.length;
};
/** The longest match from a start by a JavaScript `RegExp` tried at each length (the test's oracle only). */
const regexAt = (pattern: string, max: number) => {
  const re = new RegExp(`^(?:${pattern})$`, 'i');
  return (text: string, s: number) => {
    for (let len = Math.min(max, text.length - s); len > 0; len--) if (re.test(text.slice(s, s + len))) return len;
    return 0;
  };
};

describe('the built-in matchers against the naive oracle', () => {
  const rand = seeded(7);
  for (let round = 0; round < 6; round++) {
    const letters = randomLetters(rand, 400 + round * 37, round % 2 ? 'ACGT' : 'ACGTN');
    const exact = randomLetters(rand, 3 + (round % 3), 'ACGT');
    const degenerate = randomLetters(rand, 4, 'ACGTRYN');
    for (const cyclic of [false, true]) {
      it(`exact ${exact} on ${letters.length} letters${cyclic ? ', a circle' : ''}`, async () => {
        const r = await runSearch(ExactMatcher, exact, [textTarget(letters, { cyclic })], {}, { chunkLetters: 64 });
        expect(r.hits).toEqual(oracle(letters, exactAt(exact), { cyclic }));
      });
      it(`degenerate ${degenerate} on ${letters.length} letters${cyclic ? ', a circle' : ''}`, async () => {
        const r = await runSearch(
          DegenerateMatcher,
          degenerate,
          [textTarget(letters, { cyclic })],
          {},
          { chunkLetters: 50 },
        );
        expect(r.hits).toEqual(oracle(letters, degenerateAt(degenerate), { cyclic }));
      });
    }
    for (const pattern of ['GG.{2,8}CC', 'A(C|GT)+A', 'T[AG]{2}C?', 'C(AA|A)*G']) {
      it(`regex ${pattern} on ${letters.length} letters`, async () => {
        const r = await runSearch(RegexMatcher, pattern, [textTarget(letters)], {}, { chunkLetters: 90 });
        expect(r.hits).toEqual(oracle(letters, regexAt(pattern, 1000)));
      });
    }
  }
  it('a long pattern (multiword shift-and, 70 positions) finds what the oracle does', async () => {
    const rand2 = seeded(3);
    const unit = randomLetters(rand2, 70, 'ACGT');
    const letters = randomLetters(rand2, 200, 'ACGT') + unit + randomLetters(rand2, 50, 'ACGT') + unit;
    const r = await runSearch(DegenerateMatcher, unit.slice(0, 33) + 'N' + unit.slice(34), [textTarget(letters)], {
      strands: 'forward',
    });
    expect(r.hits.map((h) => h[1])).toEqual([200, 320]);
  });
});

describe('the oracle at small chunks: no hit lost or added at a boundary, either strand', () => {
  const rand = seeded(11);
  // long runs of A and T, so matches cross the chunks on both strands
  const texts = [0, 1, 2].map((k) => randomLetters(rand, 120 + k * 13, 'AAAAATTTCG'));
  for (const chunkLetters of [5, 7, 11]) {
    for (const [k, letters] of texts.entries()) {
      for (const pattern of ['A+', 'A+C?', 'GA*', 'T(A|T)*G', 'C(AA|A)*G']) {
        it(`regex ${pattern}, text ${k}, chunks of ${chunkLetters}`, async () => {
          const r = await runSearch(RegexMatcher, pattern, [textTarget(letters)], {}, { chunkLetters });
          expect(r.hits).toEqual(oracle(letters, regexAt(pattern, 1000)));
        });
      }
      for (const cyclic of [false, true]) {
        it(`exact AAT and degenerate WWAN, text ${k}, chunks of ${chunkLetters}${cyclic ? ', a circle' : ''}`, async () => {
          const target = [textTarget(letters, { cyclic })];
          const e = await runSearch(ExactMatcher, 'AAT', target, {}, { chunkLetters });
          expect(e.hits).toEqual(oracle(letters, exactAt('AAT'), { cyclic }));
          const d = await runSearch(DegenerateMatcher, 'WWAN', target, {}, { chunkLetters });
          expect(d.hits).toEqual(oracle(letters, degenerateAt('WWAN'), { cyclic }));
        });
      }
    }
  }
});

describe('the hit rules', () => {
  it('a site the same on both strands is one hit, strand 0 (GAATTC)', async () => {
    const r = await runSearch(ExactMatcher, 'GAATTC', [textTarget('TTGAATTCAA')]);
    expect(r.hits).toEqual([[0, 2, 8, 0]]);
  });
  it('a minus-strand hit is mapped back to forward positions', async () => {
    // GGATCC's reverse complement is itself; CCA's is TGG
    const r = await runSearch(ExactMatcher, 'CCA', [textTarget('AATGGAA')]);
    expect(r.hits).toEqual([[0, 2, 5, -1]]);
    const f = await runSearch(ExactMatcher, 'CCA', [textTarget('AATGGAA')], { strands: 'forward' });
    expect(f.hits).toEqual([]);
  });
  it('EC-10, EC-9: a hit across a circle’s origin is one hit, held unrolled (end past the length)', async () => {
    // GAATTC split at the origin of a 10-letter circle: TTC at the start, GAA at the end
    const r = await runSearch(ExactMatcher, 'GAATTC', [textTarget('TTCAAAAGAA', { cyclic: true })]);
    expect(r.hits).toEqual([[0, 7, 13, 0]]);
    // the same letters, not cyclic: no hit
    const linear = await runSearch(ExactMatcher, 'GAATTC', [textTarget('TTCAAAAGAA')]);
    expect(linear.hits).toEqual([]);
  });
  it('a pattern longer than the circle has no hit on it (at most one lap)', async () => {
    const r = await runSearch(ExactMatcher, 'ACGTA', [textTarget('ACGT', { cyclic: true })], { strands: 'forward' });
    expect(r.hits).toEqual([]);
  });
  it('overlapping hits are all kept (AA in AAA is two)', async () => {
    const r = await runSearch(ExactMatcher, 'AA', [textTarget('AAA')], { strands: 'forward' });
    expect(r.hits).toEqual([
      [0, 0, 2, 1],
      [0, 1, 3, 1],
    ]);
  });
  it('case-insensitive, U as T on nucleotides', async () => {
    const r = await runSearch(ExactMatcher, 'acu', [textTarget('ggACTgg', { alphabet: 'rna' })], {
      strands: 'forward',
    });
    expect(r.hits).toEqual([[0, 2, 5, 1]]);
  });
  it('degenerate is subset matching: N in the sequence matches only N; R matches R, never A alone', async () => {
    expect(
      (
        await runSearch(DegenerateMatcher, 'GANTC', [textTarget('GAATC GANTC'.replace(' ', 'T'))], {
          strands: 'forward',
        })
      ).hits,
    ).toEqual([
      [0, 0, 5, 1],
      [0, 6, 11, 1],
    ]);
    // a sequence N under a pattern A: not a hit
    expect((await runSearch(DegenerateMatcher, 'GAATC', [textTarget('GANTC')], { strands: 'forward' })).hits).toEqual(
      [],
    );
    // a sequence R under a pattern D (A, G or T): a hit; under a pattern A: not
    expect((await runSearch(DegenerateMatcher, 'CDC', [textTarget('CRC')], { strands: 'forward' })).hits).toEqual([
      [0, 0, 3, 1],
    ]);
    expect((await runSearch(DegenerateMatcher, 'CAC', [textTarget('CRC')], { strands: 'forward' })).hits).toEqual([]);
  });
  it('protein degenerate letters: X any residue, B D or N, Z E or Q, J I or L; one strand', async () => {
    const t = textTarget('MNDSTQEIL', { alphabet: 'protein' });
    expect((await runSearch(DegenerateMatcher, 'BX', [t])).hits).toEqual([
      [0, 1, 3, 0],
      [0, 2, 4, 0],
    ]);
    expect((await runSearch(DegenerateMatcher, 'ZJ', [t])).hits).toEqual([[0, 6, 8, 0]]);
  });
  it('gaps are not letters: a row is searched over its residues, a hit spanning the gap columns between them', async () => {
    const r = await runSearch(ExactMatcher, 'GAAT', [textTarget('CG-A.AT-C', { alphabet: 'dna' })], {
      strands: 'forward',
    });
    expect(r.hits).toEqual([[0, 1, 7, 1]]);
  });
  it('EC-37: the system’s alphabet decides; a target that cannot hold the pattern is skipped, said, the others searched', async () => {
    const r = await runSearch(ExactMatcher, 'EEK', [
      textTarget('AAEEKAA', { system: '/data/dna', alphabet: 'dna' }),
      textTarget('MEEKL', { system: '/data/p', alphabet: 'protein' }),
    ]);
    expect(r.hits).toEqual([[1, 1, 4, 0]]);
    expect(r.issues.map((i) => [i.code, i.address])).toEqual([['seq.search.alphabet', '/data/dna']]);
    // U in a protein pattern under `degenerate`
    const u = await runSearch(DegenerateMatcher, 'MU', [textTarget('MU', { alphabet: 'protein' })]);
    expect(u.issues.map((i) => i.code)).toEqual(['seq.search.alphabet']);
  });
  it('hits are ordered by target, then start, end and strand; no hit spans two targets', async () => {
    const r = await runSearch(ExactMatcher, 'ACG', [textTarget('CGTACG'), textTarget('ACGAC')]);
    expect(r.hits).toEqual([
      [0, 0, 3, -1],
      [0, 3, 6, 1],
      [1, 0, 3, 1],
    ]);
  });
  it('batches are a prefix of the final order: one a chunk, each sorted, a hit in the chunk it starts in', async () => {
    const letters = 'ACGTACGTAC'.repeat(30);
    const r = await runSearch(ExactMatcher, 'GTAC', [textTarget(letters)], {}, { chunkLetters: 16 });
    expect(r.batches.length).toBe(Math.ceil(letters.length / 16));
    const starts = r.hits.map((h) => h[1]);
    expect(starts).toEqual([...starts].sort((a, b) => a - b));
    expect(r.hits.length).toBe(oracle(letters, exactAt('GTAC')).length);
  });
  it('within: only the ranges named are searched, a hit inside them', async () => {
    const r = await runSearch(ExactMatcher, 'AC', [textTarget('ACACACAC', { within: [2, 5] })], { strands: 'forward' });
    expect(r.hits).toEqual([[0, 2, 4, 1]]);
  });
  it('regex hits: one a start, the longest from it', async () => {
    const r = await runSearch(RegexMatcher, 'A+', [textTarget('CAAAC', { alphabet: 'protein' })]);
    expect(r.hits).toEqual([
      [0, 1, 4, 0],
      [0, 2, 4, 0],
      [0, 3, 4, 0],
    ]);
  });
  it('regex anchors: ^ and $ at a system’s ends; a cyclic system has none', async () => {
    const t = textTarget('ACGTAC', { alphabet: 'protein' });
    expect((await runSearch(RegexMatcher, '^AC', [t])).hits).toEqual([[0, 0, 2, 0]]);
    expect((await runSearch(RegexMatcher, 'AC$', [t])).hits).toEqual([[0, 4, 6, 0]]);
    const c = await runSearch(RegexMatcher, '^AC', [textTarget('ACGT', { cyclic: true })]);
    expect(c.issues.map((i) => i.code)).toEqual(['seq.search.pattern.anchor']);
  });
  it('a regex hit across a circle’s origin', async () => {
    const r = await runSearch(RegexMatcher, 'TTC+A', [textTarget('CCAGGGTT', { cyclic: true })], {
      strands: 'forward',
    });
    expect(r.hits).toEqual([[0, 6, 11, 1]]);
  });
});
