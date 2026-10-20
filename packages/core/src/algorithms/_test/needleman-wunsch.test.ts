/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { AlignmentOp, type GlobalAlignScoring, globalAlign } from '../needleman-wunsch.js';

const codes = (s: string) => Array.from(s, (c) => c.charCodeAt(0));
const show = (ops: ArrayLike<number>) => Array.from(ops, (o) => 'MID'[o]).join('');

const base: GlobalAlignScoring = {
  match: 5,
  mismatch: -4,
  gapOpen: 10,
  gapExtend: 1,
  leadingGaps: 'penalized',
  trailingGaps: 'penalized',
  gapTies: 'open',
};

/** The score of `ops` under `s`, walked column by column: an independent reading of the scoring rules. */
function scoreOf(a: number[], b: number[], ops: ArrayLike<number>, s: GlobalAlignScoring): number {
  const pair = (x: number, y: number) =>
    s.matrix ? s.matrix.scores[x * s.matrix.size + y] : x === y ? (s.match ?? 0) : (s.mismatch ?? 0);
  let score = 0;
  let i = 0;
  let j = 0;
  let prev = -1;
  for (let k = 0; k < ops.length; k++) {
    const op = ops[k];
    if (op === AlignmentOp.Pair) {
      score += pair(a[i++], b[j++]);
    } else {
      const leading = op === AlignmentOp.OnlyA ? j === 0 : i === 0;
      const trailing = op === AlignmentOp.OnlyA ? j === b.length : i === a.length;
      const first = prev !== op;
      if (leading) {
        if (s.leadingGaps === 'penalized') score -= first ? s.gapOpen : s.gapExtend;
        else if (s.leadingGaps === 'flat' && first) score -= s.gapOpen;
      } else if (!(trailing && s.trailingGaps === 'free')) {
        score -= first ? s.gapOpen : s.gapExtend;
      }
      if (op === AlignmentOp.OnlyA) i++;
      else j++;
    }
    prev = op;
  }
  expect([i, j]).toEqual([a.length, b.length]);
  return score;
}

/** Every alignment of `a` and `b`, as op arrays. */
function* allAlignments(n: number, m: number, prefix: number[] = []): Generator<number[]> {
  if (n === 0 && m === 0) yield prefix;
  if (n > 0 && m > 0) yield* allAlignments(n - 1, m - 1, [...prefix, AlignmentOp.Pair]);
  if (n > 0) yield* allAlignments(n - 1, m, [...prefix, AlignmentOp.OnlyA]);
  if (m > 0) yield* allAlignments(n, m - 1, [...prefix, AlignmentOp.OnlyB]);
}

function generator(seed: number): () => number {
  let s = seed;
  return () => {
    s = (Math.imul(s, 1103515245) + 12345) >>> 0;
    return s >>> 16;
  };
}

describe('globalAlign', () => {
  it('aligns with affine gaps', () => {
    const r = globalAlign(codes('TTTTTTGGGGG'), codes('TTTGGGG'), base);
    expect(show(r.ops)).toBe('MMMIIIIMMMM');
    expect(r.score).toBe(7 * 5 - 10 - 3);
    expect(
      show(
        globalAlign(codes('TTTTTTGGGGG'), codes('TTTGGGG'), { ...base, leadingGaps: 'free', trailingGaps: 'free' }).ops,
      ),
    ).toBe('IIIMMMMMMMI');
  });

  it('covers empty sequences', () => {
    expect(show(globalAlign([], [], base).ops)).toBe('');
    expect(globalAlign([], [], base).score).toBe(0);
    expect(show(globalAlign(codes('AC'), [], base).ops)).toBe('II');
    expect(globalAlign(codes('AC'), [], base).score).toBe(-11);
    expect(globalAlign([], codes('ACG'), { ...base, leadingGaps: 'flat' }).score).toBe(-10);
  });

  it('reads a substitution matrix by code', () => {
    // codes 0 and 1; 0 pairs with 1 better than with itself
    const matrix = { size: 2, scores: [1, 8, 8, 1] };
    const r = globalAlign([0, 0], [1, 1], { ...base, matrix });
    expect(show(r.ops)).toBe('MM');
    expect(r.score).toBe(16);
  });

  it('breaks gap ties as gapTies says', () => {
    // Linear gaps: 'A' pairs with either 'A' of 'AAC'; the gaps join (extend) or open apart at equal cost.
    const linear = { ...base, gapOpen: 2, gapExtend: 2 };
    expect(show(globalAlign(codes('A'), codes('AAC'), linear).ops)).toBe('DMD');
    expect(show(globalAlign(codes('A'), codes('AAC'), { ...linear, gapTies: 'extend' }).ops)).toBe('MDD');
  });

  it('is optimal and reports its own score: every alignment of small seeded pairs', () => {
    const next = generator(7);
    const variants: GlobalAlignScoring[] = [
      base,
      { ...base, leadingGaps: 'flat', gapTies: 'extend' },
      { ...base, leadingGaps: 'free', trailingGaps: 'free' },
      { ...base, gapOpen: 3, gapExtend: 3, match: 2, mismatch: -1 },
      { ...base, gapOpen: 1.5, gapExtend: 0.25, matrix: { size: 3, scores: [3, -1, 0, -1, 2, -2, 0, -2, 4] } },
    ];
    for (let k = 0; k < 300; k++) {
      const s = variants[k % variants.length];
      const letters = s.matrix ? 3 : 2;
      const a = Array.from({ length: next() % 6 }, () => next() % letters);
      const b = Array.from({ length: next() % 6 }, () => next() % letters);
      const r = globalAlign(a, b, s);
      let best = -Infinity;
      for (const ops of allAlignments(a.length, b.length)) best = Math.max(best, scoreOf(a, b, ops, s));
      expect(r.score, `${a} | ${b} | ${k % variants.length}`).toBe(best);
      expect(scoreOf(a, b, r.ops, s)).toBe(r.score);
    }
  });
});
