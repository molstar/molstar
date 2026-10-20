/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Column } from '@molstar/core/data/db';
import { IndexPairBonds } from '../property/bonds/index-pair.js';

/** The graph's slots and edge props as plain arrays. */
function read(bonds: IndexPairBonds) {
  const { offset, a, b, edgeProps } = bonds.bonds;
  const props: Record<string, number[]> = {};
  for (const k of Object.keys(edgeProps).sort()) props[k] = Array.from((edgeProps as any)[k] as ArrayLike<number>);
  return { offset: Array.from(offset), a: Array.from(a), b: Array.from(b), props, hasOperators: bonds.hasOperators };
}

describe('IndexPairBonds', () => {
  it('stores every prop of each bond in both directions', () => {
    const bonds = IndexPairBonds.fromData({
      pairs: {
        indexA: Column.ofIntArray([0, 1, 3]),
        indexB: Column.ofIntArray([1, 2, 0]),
        key: Column.ofIntArray([10, 11, 12]),
        operatorA: Column.ofIntArray([0, 1, 2]),
        operatorB: Column.ofIntArray([3, 1, 0]),
        order: Column.ofIntArray([1, 2, 3]),
        distance: Column.ofFloatArray([1.53, -1, 2.0000001]),
        flag: Column.ofIntArray([1, 0x10 | 1, 0x20]),
      },
      count: 4,
    });
    expect(read(bonds)).toEqual({
      offset: [0, 2, 4, 5, 6],
      a: [0, 0, 1, 1, 2, 3],
      b: [1, 3, 0, 2, 1, 0],
      props: {
        distance: [1.53, 2.0000001, 1.53, -1, -1, 2.0000001],
        flag: [1, 32, 1, 17, 17, 32],
        key: [10, 12, 10, 11, 11, 12],
        operatorA: [0, 0, 3, 1, 1, 2],
        operatorB: [3, 2, 0, 1, 1, 0],
        order: [1, 3, 1, 2, 2, 3],
      },
      hasOperators: true,
    });
    expect([...bonds.bySameOperator]).toEqual([[1, [3, 4]]]);
  });

  it('fills defaults for props not given', () => {
    const bonds = IndexPairBonds.fromData({
      pairs: { indexA: Column.ofIntArray([0, 2]), indexB: Column.ofIntArray([1, 1]) },
      count: 3,
    });
    expect(read(bonds)).toEqual({
      offset: [0, 1, 3, 4],
      a: [0, 1, 1, 2],
      b: [1, 0, 2, 1],
      props: {
        distance: [-1, -1, -1, -1],
        flag: [1, 1, 1, 1],
        key: [-1, -1, -1, -1],
        operatorA: [-1, -1, -1, -1],
        operatorB: [-1, -1, -1, -1],
        order: [1, 1, 1, 1],
      },
      hasOperators: false,
    });
  });
});
