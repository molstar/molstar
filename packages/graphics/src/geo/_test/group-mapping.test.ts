/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createGroupMapping } from '../util.js';

function groupsOf(mapping: ReturnType<typeof createGroupMapping>, group: number) {
  const { indices, offsets } = mapping;
  return Array.from(indices).slice(offsets[group], offsets[group + 1]);
}

describe('createGroupMapping', () => {
  it('maps every group id, also ids at or above the data count (perf-audit-2 08/B2)', () => {
    // 4 items in groups 5, 7, 5, 9: ids beyond the item count, as a shape with sparse group ids has
    const mapping = createGroupMapping([5, 7, 5, 9], 4);
    expect(groupsOf(mapping, 5)).toEqual([0, 2]);
    expect(groupsOf(mapping, 7)).toEqual([1]);
    expect(groupsOf(mapping, 9)).toEqual([3]);
    expect(groupsOf(mapping, 0)).toEqual([]);
    expect(groupsOf(mapping, 6)).toEqual([]);
  });

  it('reads one group per step', () => {
    // 2 items of 3 entries each (as cylinders repeat a group per vertex)
    const mapping = createGroupMapping([0, 0, 0, 1, 1, 1], 2, 3);
    expect(groupsOf(mapping, 0)).toEqual([0]);
    expect(groupsOf(mapping, 1)).toEqual([3]);
  });
});
