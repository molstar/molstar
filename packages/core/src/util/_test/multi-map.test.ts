/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { MultiMap } from '../multi-map.js';

describe('MultiMap', () => {
  it('collects the values of a key in insertion order, missing keys have none', () => {
    const map = new MultiMap<string, number>();
    expect(map.get('a')).toBeUndefined();
    map.add('a', 1);
    map.add('b', 2);
    map.add('a', 3);
    expect(map.get('a')).toEqual([1, 3]);
    expect(map.get('b')).toEqual([2]);
    expect(map.get('c')).toBeUndefined();
  });
});
