/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What the axis reads of a path (_spec/README.md "Resolve", Paths): points place by place (each nucleotide of a split
// codon its own edge), a window's end, the majority direction, the coding runs; on every encoding of the same map.

import { describe, expect, it } from 'vitest';
import { codonLinear } from '../codon.js';
import { expand, type Linear, RunsLinear, toRuns } from '../linear.js';
import {
  codingRanges,
  composeEdges,
  inversePath,
  type Path,
  pathDirection,
  projectEnd,
  projectPoint,
} from '../paths.js';

const codonPath = (l: Linear): Path =>
  composeEdges('/p', [{ map: '/m', from: '/p', to: '/g', backward: false, step: 1, linear: l }])!;
const stepPath = (l: Linear): Path =>
  composeEdges('/a', [{ map: '/m', from: '/a', to: '/b', backward: false, step: 0, linear: l }])!;

// exons [10, 14) and [20, 25) on the minus strand: 9 nucleotides, 3 residues; residue 0 is 24, 23, 22; residue 1 is
// 21, 20 and then 13 (split by the intron); residue 2 is 12, 11, 10
const minus = codonLinear([10, 20], [14, 25], '-');
const plus = codonLinear([10, 20], [14, 25], '+');
const all = (l: Linear) => [l, toRuns(l), expand(l, 0, 9)].map(codonPath);

describe('points through a codon step', () => {
  it('each nucleotide of a split codon is an edge of its own', () => {
    for (const p of all(minus)) {
      // residue 1 starts at 21's far edge (22), runs down to 20, then jumps the intron to 13's far edge (14)
      expect(projectPoint(p, 1)).toBe(22);
      expect(projectPoint(p, 3.5 / 3)).toBeCloseTo(21.5);
      expect(projectPoint(p, 4.5 / 3)).toBeCloseTo(20.5);
      expect(projectPoint(p, 5.001 / 3)).toBeCloseTo(14, 2);
      expect(projectPoint(p, 5.5 / 3)).toBeCloseTo(13.5);
      // every nucleotide of the coding sequence is some place's image
      const reached: number[] = [];
      for (let f = 0; f < 9; f++) reached.push(Math.floor(projectPoint(p, (f + 0.5) / 3)!));
      expect(reached.sort((a, b) => a - b)).toEqual([10, 11, 12, 13, 20, 21, 22, 23, 24]);
    }
    for (const p of all(plus)) {
      expect(projectPoint(p, 4.5 / 3)).toBeCloseTo(20.5);
      expect(projectPoint(p, 3.5 / 3)).toBeCloseTo(13.5);
    }
  });

  it('a residue a hair below a whole third is that third: never the place before it, across the intron', () => {
    // a hair below 5 / 3: three times it is 4.999…, not 5 (a window's edge read back through offsets)
    const x = 5 / 3 - 1e-15;
    expect(3 * x).toBeLessThan(5);
    for (const p of all(minus)) {
      // place 5 is nucleotide 13 (its edge 14), not the end of place 4 (20, across the intron)
      expect(projectPoint(p, x)).toBeCloseTo(14, 9);
      expect(projectPoint(p, 5 / 3)).toBe(14);
    }
  });

  it('a window end is the far edge of its last place', () => {
    for (const p of all(minus)) {
      expect(projectEnd(p, 2)).toBe(13);
      expect(projectEnd(p, 1)).toBe(22);
      expect(projectEnd(p, 1.5)).toBeCloseTo(projectPoint(p, 1.5)!);
      // back: a nucleotide window's end onto residues, a third a nucleotide
      const back = inversePath(p);
      expect(projectEnd(back, 22)).toBeCloseTo(1);
      expect(projectPoint(back, 24.5)).toBeCloseTo(0.5 / 3);
    }
  });

  it('the coding runs are the nucleotide side, either way', () => {
    for (const p of all(minus)) {
      expect(codingRanges(p)).toEqual([10, 14, 20, 25]);
      expect(codingRanges(inversePath(p))).toEqual([10, 14, 20, 25]);
    }
    expect(codingRanges(stepPath(new RunsLinear([{ from: 0, to: 0, length: 5, dir: 1 }])))).toBeUndefined();
  });
});

describe('a path’s direction', () => {
  it('is the majority of its runs’ lengths, mixed when they go both ways', () => {
    for (const p of all(minus)) expect(pathDirection(p)).toEqual({ dir: -1, mixed: false });
    // an inversion of 3 inside 10 forward
    const inversion = new RunsLinear([
      { from: 0, to: 0, length: 4, dir: 1 },
      { from: 4, to: 6, length: 3, dir: -1 },
      { from: 7, to: 7, length: 3, dir: 1 },
    ]);
    expect(pathDirection(stepPath(inversion))).toEqual({ dir: 1, mixed: true });
    // mostly reversed, its first run forward: the majority, not the first run
    const mostly = new RunsLinear([
      { from: 0, to: 0, length: 2, dir: 1 },
      { from: 2, to: 20, length: 10, dir: -1 },
    ]);
    expect(pathDirection(stepPath(mostly))).toEqual({ dir: -1, mixed: true });
  });
});
