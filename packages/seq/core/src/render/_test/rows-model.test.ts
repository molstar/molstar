/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { matchLeaves, readNewick } from '@molstar/seqx-core/data/alignment-readers';
import { describe, expect, it } from 'vitest';
import {
  allClades,
  buildClades,
  cladesAtDepth,
  cladesHolding,
  cladeNamed,
  cladesWithin,
  columnMeasures,
  consensus,
  DisplayRows,
  type HeldCells,
  logoHeights,
  measureEntries,
  profileOf,
} from '../rows-model.js';

function held(letters: string[], labels = letters.map((_, i) => `r${i}`)): HeldCells {
  const C = letters[0].length;
  const codes = new Uint8Array(letters.length * C);
  letters.forEach((l, i) => {
    for (let j = 0; j < C; j++) codes[i * C + j] = l.charCodeAt(j);
  });
  return { rows: letters.length, columns: C, codes, labels, first: new Int32Array(letters.length).fill(1), version: 1 };
}

/** A random binary tree over `n` leaves named r0..r(n−1) in a shuffled order, with lengths. */
function randomTree(rng: PCG, n: number): string {
  const names = Array.from({ length: n }, (_, i) => `r${i}`);
  for (let i = n - 1; i > 0; i--) {
    const j = rng.int() % (i + 1);
    [names[i], names[j]] = [names[j], names[i]];
  }
  const build = (a: number, b: number): string => {
    if (b - a === 1) return `${names[a]}:${(rng.int() % 10) / 10}`;
    const m = a + 1 + (rng.int() % (b - a - 1));
    return `(${build(a, m)},${build(m, b)}):${(rng.int() % 10) / 10}`;
  };
  return `${build(0, n)};`;
}

describe('the rows’ model', () => {
  const labels = ['a', 'b', 'c', 'd', 'e'];
  const tree = readNewick('((a:1,b:1)ab:1,((c:1,d:2):1,e:1):1)root;');
  const clades = buildClades(tree, matchLeaves(tree, labels.map((label) => ({ label }))), 5)!;

  it('clades: named by their label, else first..last; runs of leaves; heights by branch length', () => {
    expect(clades.name[0]).toBe('root');
    const cd = cladeNamed(clades, 'c..d');
    expect(cd).toBeGreaterThan(0);
    expect([clades.lo[cd], clades.hi[cd]]).toEqual([2, 4]);
    expect(clades.height[cd]).toBe(2);
    expect(cladesAtDepth(clades, 1).map((c) => clades.name[c])).toEqual(['ab', 'c..e']);
    expect(cladesWithin(clades, 1).map((c) => clades.name[c])).toEqual(['ab']);
    expect(cladesWithin(clades, 2).map((c) => clades.name[c])).toEqual(['ab', 'c..d']);
    expect(allClades(clades).map((c) => clades.name[c])).toEqual(['ab', 'c..e', 'c..d']);
    expect(cladesHolding(clades, 3).map((c) => clades.name[c])).toEqual(['root', 'c..e', 'c..d']);
  });

  it('display rows: a folded clade one row, a pinned row in the band and in place, hidden rows out and counted', () => {
    const d = new DisplayRows();
    const cd = cladeNamed(clades, 'c..d');
    d.set(5, undefined, clades, [cd], [4], new Set([1]), true);
    // the pinned band: e; then a, (c..d), e
    expect(d.pinned).toBe(1);
    expect(d.count).toBe(4);
    expect(Array.from(d.entry.subarray(0, 4))).toEqual([4, 0, -1 - cd, 4]);
    expect(Array.from(d.members.get(2)!)).toEqual([2, 3]);
    expect(d.displayRow(3)).toBe(2);
    expect(d.displayRow(1)).toBe(-1);
    expect(d.displayRow(4)).toBe(0);
    expect(d.hidden).toBe(1);
    // a clade with one shown row folds to that row
    d.set(5, undefined, clades, [cladeNamed(clades, 'ab')], [], new Set([0]), true);
    expect(Array.from(d.entry.subarray(0, d.count))).toEqual([1, 2, 3, 4]);
    d.setIdentity(1_000_000);
    expect(d.rowAt(123_456)).toBe(123_456);
    expect(d.displayRow(9)).toBe(9);
  });

  it('property: over random trees and folds, every shown row is placed once and found where it is drawn', () => {
    const rng = new PCG(21);
    const d = new DisplayRows();
    for (let t = 0; t < 60; t++) {
      const n = 2 + (rng.int() % 40);
      const tr = readNewick(randomTree(rng, n));
      const c = buildClades(tr, matchLeaves(tr, Array.from({ length: n }, (_, i) => ({ label: `r${i}` }))), n)!;
      const folds = allClades(c).filter(() => rng.float() < 0.3);
      const hidden = new Set(Array.from({ length: n }, (_, i) => i).filter(() => rng.float() < 0.1));
      d.set(n, undefined, c, folds, [], hidden, true);
      const seen = new Map<number, number>();
      for (let i = 0; i < d.count; i++) {
        const e = d.rowAt(i);
        const rows = e >= 0 ? [e] : Array.from(d.members.get(i)!);
        for (const r of rows) {
          expect(seen.has(r)).toBe(false);
          seen.set(r, i);
        }
      }
      for (let r = 0; r < n; r++) {
        if (hidden.has(r)) expect(d.displayRow(r)).toBe(-1);
        else expect(d.displayRow(r)).toBe(seen.get(r));
      }
      expect(seen.size + hidden.size).toBe(n);
    }
  });

  it('measures, consensus and the profile from the held cells', () => {
    const h = held(['ACD-', 'ACE-', 'GCDA', 'ACDA']);
    const m = columnMeasures(h, 0);
    expect(Array.from(m.share)).toEqual([0.75, 1, 0.75, 1]);
    expect(m.identity[0]).toBe(0.75);
    expect(Number.isNaN(m.identity[3])).toBe(true);
    expect(m.rowsValue[1]).toBeCloseTo(2 / 3);
    const entries = measureEntries(h, m);
    const datum = { field: () => undefined, value: 0, letter: 'G', start: 0, end: 1 };
    expect(entries.get('matches-reference')!.read(datum, {})).toBe(0);
    expect(entries.get('share')!.read(datum, {})).toBe(0.75);
    const out = consensus(h, [0, 1, 2], new Uint8Array(12));
    expect(String.fromCharCode(out[0], out[3], out[6], out[9])).toBe('ACD-');
    expect(out[1]).toBe(170);
    const p = profileOf(h);
    expect(String.fromCharCode(...p.letters)).toBe('ACDEG');
    const hs = logoHeights(p, 1, 'bits', 'none', 20, new Float32Array(5));
    expect(hs[1]).toBeCloseTo(Math.log2(20));
    const pr = logoHeights(p, 0, 'probability', 'none', 20, new Float32Array(5));
    expect(pr[0]).toBeCloseTo(0.75);
    // the small-sample correction lowers it
    expect(logoHeights(p, 1, 'bits', 'small-sample', 20, new Float32Array(5))[1]).toBeLessThan(hs[1]);
  });
});
