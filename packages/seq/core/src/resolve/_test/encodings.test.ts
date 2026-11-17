/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The edge cases of maps on every encoding of the same map (implicit where it is arithmetic, runs, expanded over its
// domain), with equal answers; and seeded properties at small counts.

import { PCG } from '@molstar/core/data/util/hash-functions';
import { describe, expect, it } from 'vitest';
import { codonLinear, codonPairs, codonRuns, codonStatus, translateCodon } from '../codon.js';
import {
  canonicalRuns,
  composeLinear,
  expand,
  ImplicitLinear,
  type Linear,
  linearsEqual,
  type Run,
  RunsLinear,
  toRuns,
} from '../linear.js';
import { composeEdges, type Edge, type Path, pathStatus, projectBoundary, projectRanges } from '../paths.js';

/** The same map in each encoding it has: the given one, runs, and expanded over `[lo, hi)`. */
function encodings(l: Linear, lo: number, hi: number): Linear[] {
  return [l, toRuns(l), expand(l, lo, hi)];
}

/** Every encoding of a codon map as a one-step path from the protein (`p`) to the gene (`g`). */
function codonPaths(linear: Linear, codingLength: number): Path[] {
  return encodings(linear, 0, codingLength).map((l) => {
    const edge: Edge = { map: '/data/m', from: '/p', to: '/g', backward: false, step: 1, linear: l };
    return composeEdges('/p', [edge])!;
  });
}

function backPath(p: Path): Path {
  const e = p.edges[0];
  return composeEdges('/g', [{ ...e, from: '/g', to: '/p', backward: true, step: -1, linear: e.linear.inverse() }])!;
}

/** Image ranges as sorted, merged pairs. */
function merged(pairs: number[]): number[] {
  const order: [number, number][] = [];
  for (let i = 0; i < pairs.length; i += 2) order.push([pairs[i], pairs[i + 1]]);
  order.sort((a, b) => a[0] - b[0]);
  const out: number[] = [];
  for (const [s, e] of order) {
    if (out.length && s <= out[out.length - 1]) out[out.length - 1] = Math.max(out[out.length - 1], e);
    else out.push(s, e);
  }
  return out;
}

const same = <T>(values: T[]) => {
  for (const v of values) expect(v).toEqual(values[0]);
  return values[0];
};

// the p53 example's gene: exons [100, 150) and [160, 200), 90 nucleotides, 30 residues
const Starts = [100, 160];
const Ends = [150, 200];

describe('maps: offset on a cyclic target (EC-9 to EC-11)', () => {
  // six positions onto a plasmid of 10 at +7: 0..2 → 7..9, 3..5 → 0..2
  const offset = new ImplicitLinear(7, 1, 0, 6, 10);

  it('a segment crossing the origin is split when built (EC-11)', () => {
    expect(offset.runs()).toEqual([
      { from: 0, to: 7, length: 3, dir: 1 },
      { from: 3, to: 0, length: 3, dir: 1 },
    ]);
  });

  it('a wrapped image is two runs, never one past the end, on every encoding (EC-9, EC-10)', () => {
    const fwd = same(encodings(offset, 0, 6).map((l) => merged(runsToPairs(l.forward(1, 5)))));
    expect(fwd).toEqual([0, 2, 8, 10]);
    const back = same(
      encodings(offset, 0, 6).map((l) => merged([...runsToPairs(l.backward(8, 10)), ...runsToPairs(l.backward(0, 1))])),
    );
    expect(back).toEqual([1, 4]);
  });

  it('what lies past a non-cyclic target maps nowhere, and is no error (EC-21)', () => {
    const partial = new ImplicitLinear(-5, 1, 5, 10);
    const answers = encodings(partial, 0, 10).map((l) => runsToPairs(l.forward(0, 10)));
    expect(same(answers)).toEqual([0, 5]);
    expect(same(encodings(partial, 0, 10).map((l) => runsToPairs(l.forward(0, 5))))).toEqual([]);
  });
});

function runsToPairs(runs: readonly Run[]): number[] {
  return runs.flatMap((r) => (r.dir === 1 ? [r.to, r.to + r.length] : [r.to - r.length + 1, r.to + 1]));
}

describe('maps: codon (EC-15 to EC-18)', () => {
  const minus = codonLinear(Starts, Ends, '-');
  const plus = codonLinear(Starts, Ends, '+');

  it('minus strand: the coding sequence starts at the highest exon and reads downward (EC-15)', () => {
    expect(codonRuns(Starts, Ends, '-')).toEqual([
      { from: 0, to: 199, length: 40, dir: -1 },
      { from: 40, to: 149, length: 50, dir: -1 },
    ]);
    // residue 0 is 199, 198, 197; residue 13 is split by the junction: 160, then 149 and 148
    for (const p of codonPaths(minus, 90)) {
      expect(merged(projectRanges(p, [0, 1]))).toEqual([197, 200]);
      expect(merged(projectRanges(p, [13, 14]))).toEqual([148, 150, 160, 161]);
      expect(merged(projectRanges(p, [0, 30]))).toEqual([100, 150, 160, 200]);
    }
  });

  it('plus strand: in genomic order, a junction splitting residue 16 (EC-15)', () => {
    for (const p of codonPaths(plus, 90)) {
      expect(merged(projectRanges(p, [0, 1]))).toEqual([100, 103]);
      expect(merged(projectRanges(p, [16, 17]))).toEqual([148, 150, 160, 161]);
    }
  });

  it('a residue is three nucleotides, a nucleotide one residue, a UTR or intron nucleotide none (EC-18)', () => {
    for (const p of [...codonPaths(minus, 90), ...codonPaths(plus, 90)]) {
      const back = backPath(p);
      expect(projectRanges(back, [150, 160])).toEqual([]);
      expect(projectRanges(back, [0, 100])).toEqual([]);
      expect(pathStatus(p, 0)).toBe('paired');
      expect(pathStatus(back, 155)).toBe('unpaired');
    }
    for (const p of codonPaths(minus, 90)) expect(merged(projectRanges(backPath(p), [160, 161]))).toEqual([13, 14]);
    for (const p of codonPaths(plus, 90)) expect(merged(projectRanges(backPath(p), [160, 161]))).toEqual([16, 17]);
  });

  it('a codon map within one exon is arithmetic, with the same answers as its runs', () => {
    const one = codonLinear([10], [40], '-');
    expect(one.encoding).toBe('implicit');
    for (const p of codonPaths(one, 30)) {
      expect(merged(projectRanges(p, [0, 1]))).toEqual([37, 40]);
      expect(merged(projectRanges(backPath(p), [10, 11]))).toEqual([9, 10]);
    }
  });

  it('codon status: T and U alike, ambiguity X, a first start as M, Sec and Pyl, stops pair with nothing (EC-16)', () => {
    expect(translateCodon('AUG')).toBe('M');
    expect(translateCodon('ATG')).toBe('M');
    expect(translateCodon('ANG')).toBe('X');
    expect(translateCodon('TAA')).toBe('*');
    expect(codonPairs('M', 'TTG', true)).toBe(true);
    expect(codonPairs('M', 'TTG', false)).toBe(false);
    expect(codonPairs('U', 'TGA', false)).toBe(true);
    expect(codonPairs('O', 'TAG', false)).toBe(true);
    expect(codonPairs('W', 'TAA', false)).toBe(false);
    expect(codonPairs('K', 'NNN', false)).toBe(true);
    // a gene of 12 on the plus strand: M K (conflict: TAA) W, the last residue's triplet cut short
    const l = codonLinear([0], [11], '+');
    expect(codonStatus(l, 'MKWQ', 'ATGAAATAATG')).toEqual(['paired', 'paired', 'conflict', 'unpaired']);
    // the same protein read on the minus strand: the complement, downward
    const rc = 'CATTATTTCAT';
    expect(codonStatus(codonLinear([0], [11], '-'), 'MKWQ', rc)).toEqual(['paired', 'paired', 'conflict', 'unpaired']);
  });
});

describe('maps: boundaries (EC-28, EC-29)', () => {
  it('a residue boundary lands before its residue’s first nucleotide, else after the last one before it (EC-28)', () => {
    for (const p of codonPaths(codonLinear(Starts, Ends, '+'), 90)) {
      expect(projectBoundary(p, 0)).toBe(100);
      expect(projectBoundary(p, 17)).toBe(161);
      // past the last residue: after residue 29's last nucleotide
      expect(projectBoundary(p, 30)).toBe(200);
      expect(projectBoundary(p, 40)).toBeUndefined();
    }
    for (const p of codonPaths(codonLinear(Starts, Ends, '-'), 90)) {
      expect(projectBoundary(p, 0)).toBe(200);
      expect(projectBoundary(p, 30)).toBe(100);
    }
  });

  it('a nucleotide boundary inside a codon lands after its residue (EC-29)', () => {
    for (const p of codonPaths(codonLinear(Starts, Ends, '+'), 90)) {
      const back = backPath(p);
      expect(projectBoundary(back, 100)).toBe(0);
      expect(projectBoundary(back, 101)).toBe(1);
      expect(projectBoundary(back, 103)).toBe(1);
      // at an exon's end: after the last nucleotide before it (c 49, inside residue 16); inside the intron: nowhere
      expect(projectBoundary(back, 150)).toBe(17);
      expect(projectBoundary(back, 155)).toBeUndefined();
    }
    for (const p of codonPaths(codonLinear(Starts, Ends, '-'), 90)) {
      const back = backPath(p);
      expect(projectBoundary(back, 200)).toBe(0);
      expect(projectBoundary(back, 199)).toBe(1);
      expect(projectBoundary(back, 197)).toBe(1);
    }
  });

  it('an offset boundary maps nowhere when neither neighbor maps', () => {
    for (const l of encodings(new ImplicitLinear(-5, 1, 5, 10), 0, 10)) {
      const p = composeEdges('/a', [{ map: '/m', from: '/a', to: '/b', backward: false, step: 0, linear: l }])!;
      expect(projectBoundary(p, 3)).toBeUndefined();
      expect(projectBoundary(p, 5)).toBe(0);
      expect(projectBoundary(p, 10)).toBe(5);
    }
  });
});

describe('maps: one-to-many backward (EC-20)', () => {
  // a GFP–GFP fusion of 20 onto GFP of 10: two runs onto the same positions
  const fusion = new RunsLinear([
    { from: 0, to: 0, length: 10, dir: 1 },
    { from: 10, to: 0, length: 10, dir: 1 },
  ]);

  it('reads backward as two runs, on runs and expanded alike', () => {
    for (const l of [fusion, expand(fusion, 0, 20)]) {
      expect(merged(runsToPairs(l.backward(2, 4)))).toEqual([2, 4, 12, 14]);
      expect(merged(runsToPairs(l.forward(9, 11)))).toEqual([0, 1, 9, 10]);
    }
  });
});

describe('maps: canonical runs and composition', () => {
  it('runs join position by position, whatever pieces they come in', () => {
    const a = canonicalRuns([
      { from: 0, to: 5, length: 1, dir: 1 },
      { from: 1, to: 4, length: 1, dir: 1 },
      { from: 2, to: 3, length: 2, dir: -1 },
    ]);
    expect(a).toEqual([{ from: 0, to: 5, length: 4, dir: -1 }]);
  });

  it('two implicit maps compose to an implicit one; implicit through runs to runs', () => {
    const a = new ImplicitLinear(3, 1, 0, 10);
    const b = new ImplicitLinear(20, -1, 0, 15);
    const c = composeLinear(a, b);
    expect(c.encoding).toBe('implicit');
    expect(linearsEqual(toRuns(c), composeLinear(toRuns(a), toRuns(b)))).toBe(true);
    expect(c.forward(0, 1)).toEqual([{ from: 0, to: 17, length: 1, dir: -1 }]);
    expect(composeLinear(a, toRuns(b)).encoding).toBe('runs');
  });
});

// ---- properties ----

/** A random injective map: a few runs over [0, 40) onto disjoint places of [0, 200), either direction. */
function randomRuns(rng: PCG): Run[] {
  const runs: Run[] = [];
  let from = 0;
  const used: [number, number][] = [];
  const n = 1 + (rng.int() % 4);
  for (let i = 0; i < n; i++) {
    from += rng.int() % 3;
    const length = 1 + (rng.int() % 8);
    const lo = 40 * i + (rng.int() % 20);
    used.push([lo, lo + length]);
    const dir = rng.int() % 2 ? 1 : -1;
    runs.push({ from, to: dir === 1 ? lo : lo + length - 1, length, dir });
    from += length;
  }
  return runs;
}

describe('maps: properties at small seeds', () => {
  it('forward then back is the identity where the map is total, on every encoding', () => {
    const rng = new PCG(5);
    for (let c = 0; c < 30; c++) {
      const runs = new RunsLinear(randomRuns(rng));
      for (const l of encodings(runs, 0, runs.hi)) {
        for (const r of runs.list)
          for (let k = r.from; k < r.from + r.length; k++) {
            const img = l.forward(k, k + 1);
            expect(img.length).toBe(1);
            const back = l.backward(img[0].to, img[0].to + 1);
            expect(back.map((b) => b.to)).toEqual([k]);
          }
      }
    }
  });

  it('a path equals the composition of its steps', () => {
    const rng = new PCG(11);
    for (let c = 0; c < 30; c++) {
      const a = new RunsLinear(randomRuns(rng));
      const shift = rng.int() % 50;
      const b = new ImplicitLinear(-shift, 1, shift, 250);
      const edges: Edge[] = [
        { map: '/a', from: '/x', to: '/y', backward: false, step: 0, linear: a },
        { map: '/b', from: '/y', to: '/z', backward: false, step: 0, linear: b },
      ];
      const path = composeEdges('/x', edges)!;
      for (let k = 0; k < 6; k++) {
        const s = rng.int() % 40;
        const e = s + 1 + (rng.int() % 10);
        const stepwise = runsToPairs(a.forward(s, e)).flatMap((x, i, all) =>
          i % 2 ? [] : runsToPairs(b.forward(x, all[i + 1])),
        );
        expect(merged(projectRanges(path, [s, e]))).toEqual(merged(stepwise));
      }
    }
  });
});
