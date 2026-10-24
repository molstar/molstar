/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { LRUCache } from '../lru-cache.js';

describe('LRUCache', () => {
  it('set of an existing key replaces its value (perf-audit-2 01/B4)', () => {
    const cache = LRUCache.create<number>(2);
    expect(LRUCache.set(cache, 'a', 1)).toBeUndefined();
    expect(LRUCache.set(cache, 'a', 2)).toBe(1);
    expect(cache.entries.count).toBe(1);
    expect(LRUCache.get(cache, 'a')).toBe(2);
    // the same value again is not reported as replaced
    expect(LRUCache.set(cache, 'a', 2)).toBeUndefined();
    // 'a' is the most recent: 'b' then 'c' evict 'b' first
    LRUCache.set(cache, 'b', 3);
    LRUCache.set(cache, 'a', 4);
    expect(LRUCache.set(cache, 'c', 5)).toBe(3);
    expect(LRUCache.get(cache, 'a')).toBe(4);
    expect(cache.entries.count).toBe(2);
  });
});
