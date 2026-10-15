/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { mergeParams } from '../params.js';

// The params every area patches with (params.ts).

describe('mergeParams', () => {
  it('merges a plain group one level and replaces anything else', () => {
    const base = { group: { a: 1, b: 2 }, list: [1, 2], n: 1 };
    expect(mergeParams(base, { group: { b: 3 }, list: [9], n: 2 })).toEqual({ group: { a: 1, b: 3 }, list: [9], n: 2 });
    expect(base.group).toEqual({ a: 1, b: 2 });
  });

  it('replaces a value with a kind as a whole: a relative iso value patched to absolute keeps no stale field', () => {
    const base = { alpha: 1, isoValue: { kind: 'relative', relativeValue: 1.5 } };
    const next = mergeParams(base, { isoValue: { kind: 'absolute', absoluteValue: 0.4 } });
    expect(next.isoValue).toEqual({ kind: 'absolute', absoluteValue: 0.4 });
    expect('relativeValue' in next.isoValue).toBe(false);
    expect(next.alpha).toBe(1);
  });
});
