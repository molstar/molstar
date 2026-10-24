/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { SortedArray } from '@molstar/core/data/int';
import type { ElementIndex } from '../../../../model/indexing.js';
import type { IntraUnitBonds } from '../data.js';
import { ElementSetIntraBondCache } from '../element-set-intra-bond-cache.js';

const set = (i: number) => SortedArray.ofSortedArray<ElementIndex>([i, i + 1, i + 5] as ElementIndex[]);
const bonds = (i: number) => ({ i }) as unknown as IntraUnitBonds;

describe('ElementSetIntraBondCache', () => {
  it('keeps the element sets used last only (perf-audit-2 05/P11)', () => {
    const cache = new ElementSetIntraBondCache();
    const max = ElementSetIntraBondCache.MaxEntries;
    for (let i = 0; i < max; i++) cache.set(set(i), bonds(i));
    // set 0 used again: set 1 is the least recent
    expect(cache.get(set(0))).toEqual(bonds(0));
    cache.set(set(max), bonds(max));
    expect(cache.get(set(1))).toBeUndefined();
    expect(cache.get(set(0))).toEqual(bonds(0));
    expect(cache.get(set(2))).toEqual(bonds(2));
    expect(cache.get(set(max))).toEqual(bonds(max));
    // replacing a set's bonds keeps one entry
    cache.set(set(2), bonds(-2));
    expect(cache.get(set(2))).toEqual(bonds(-2));
    for (let i = 0; i < 1000; i++) cache.set(set(1000 + i), bonds(i));
    expect((cache as unknown as { order: Set<unknown> }).order.size).toBe(max);
  });
});
