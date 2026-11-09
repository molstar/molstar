/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { createDataArray } from '@molstar/graphics/gl/renderable/util';
import { orderPositionGroup, PositionGroupStride } from '../center-group.js';

/** the items in the order of the levels of `strides`, as a pass over the items per level (coarsest first) */
function reference(count: number, strides: readonly number[]) {
  const level = (j: number) => {
    let k = 0;
    while (k < strides.length && j % strides[k] !== 0) ++k;
    return k;
  };
  const order: number[] = [];
  const offsets = [0];
  for (let i = 0; i <= strides.length; ++i) {
    for (let j = 0; j < count; ++j) if (level(j) === i) order.push(j);
    if (i < strides.length) offsets.push(order.length * 6);
  }
  return { order, offsets: strides.length ? offsets : undefined };
}

/** `count` items, the group of item j being j (the center j, 2j, 3j) */
function createItems(count: number) {
  const array = new Float32Array(count * PositionGroupStride);
  for (let j = 0; j < count; ++j) array.set([j, 2 * j, 3 * j, j], j * PositionGroupStride);
  return ValueCell.create(createDataArray(count, PositionGroupStride, Float32Array, array));
}

/** the groups in the cell's order, each item checked to be whole (its center moved with its group) */
function groups(cell: ReturnType<typeof createItems>, count: number) {
  const { array } = cell.ref.value;
  const out: number[] = [];
  for (let p = 0; p < count; ++p) {
    const g = array[p * PositionGroupStride + 3];
    expect([array[p * 4], array[p * 4 + 1], array[p * 4 + 2]]).toEqual([g, 2 * g, 3 * g]);
    out.push(g);
  }
  return out;
}

describe('orderPositionGroup', () => {
  const count = 1003;

  it('orders the items into the levels of the strides, as a pass per level', () => {
    for (const strides of [[16, 4, 1], [8, 2], [3], [1], [7, 5, 3]]) {
      const cell = createItems(count);
      const offsets = orderPositionGroup(cell, count, strides);
      const expected = reference(count, strides);
      expect(groups(cell, count)).toEqual(expected.order);
      expect(offsets).toEqual(expected.offsets);
    }
  });

  it('orders again from the order of the strides applied before, and back without strides', () => {
    const cell = createItems(count);
    orderPositionGroup(cell, count, [16, 4, 1]);
    const offsets = orderPositionGroup(cell, count, [8, 2]);
    expect(groups(cell, count)).toEqual(reference(count, [8, 2]).order);
    expect(offsets).toEqual(reference(count, [8, 2]).offsets);
    expect(orderPositionGroup(cell, count, [])).toBeUndefined();
    expect(groups(cell, count)).toEqual(reference(count, []).order);
  });

  it('takes strides that are no positive integers as before', () => {
    for (const strides of [
      [2.5, 1],
      [0, 3],
    ]) {
      const cell = createItems(count);
      const offsets = orderPositionGroup(cell, count, strides);
      expect(groups(cell, count)).toEqual(reference(count, strides).order);
      expect(offsets).toEqual(reference(count, strides).offsets);
    }
  });
});
