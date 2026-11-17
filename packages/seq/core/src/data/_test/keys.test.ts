/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { dataKey, fitLevel, nodeKey, tileRegion, tileSpan } from '../keys.js';

describe('keys', () => {
  it('are written as the protocol says', () => {
    expect(dataKey('/data/genes', 'chr1', 10, 3)).toBe('/data/genes|chr1|10|3');
    expect(nodeKey('/data/genes', 'T1')).toBe('/data/genes#T1');
    expect(nodeKey('/data/genes/(node :key T1)', 'domains')).toBe('/data/genes/(node :key T1)#domains');
  });

  it('never collide across parts holding a separator', () => {
    expect(dataKey('a|b', 'c', 1, 0)).not.toBe(dataKey('a', 'b|c', 1, 0));
    expect(nodeKey('a#b', 'c')).not.toBe(nodeKey('a', 'b#c'));
    expect(nodeKey('a\\', '#b')).not.toBe(nodeKey('a\\#', 'b'));
  });
});

describe('tiles', () => {
  it('snap a region to whole tiles of a level', () => {
    expect(tileSpan(0, 1024, 1)).toEqual([0, 0]);
    expect(tileSpan(1000, 1025, 1)).toEqual([0, 1]);
    expect(tileSpan(20_480, 30_000, 10)).toEqual([2, 2]);
    const [first, last] = tileSpan(5, 5, 1);
    expect(last).toBeLessThan(first);
  });

  it('cover positions at the level’s bin width, clipped to the length', () => {
    expect(tileRegion(2, 10)).toEqual({ start: 20_480, end: 30_720 });
    expect(tileRegion(1, 1, 1024, 1500)).toEqual({ start: 1024, end: 1500 });
    expect(tileRegion(0, 1, 4)).toEqual({ start: 0, end: 4 });
  });
});

describe('levels (EC-170)', () => {
  it('fit the coarsest listed bin no wider than a device pixel', () => {
    const levels = [1, 10, 100];
    expect(fitLevel(levels, 0.25)).toBe(1);
    expect(fitLevel(levels, 1)).toBe(1);
    expect(fitLevel(levels, 37)).toBe(10);
    expect(fitLevel(levels, 100)).toBe(100);
    expect(fitLevel(levels, 1e6)).toBe(100);
  });
});
