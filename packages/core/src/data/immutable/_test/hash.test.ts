/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { DefaultConfig, resolveConfig, sameConfig, sameValueZero } from '../config.js';
import { defaultHash, fmix32, lookupDefaultHash } from '../hash.js';
import { HashMap } from '../hash-map.js';
import { HashSet } from '../hash-set.js';

describe('defaultHash', () => {
  it('−0 hashes as 0, every NaN alike; int32 results', () => {
    expect(defaultHash(-0)).toBe(defaultHash(0));
    expect(defaultHash(NaN)).toBe(defaultHash(Number('x')));
    expect(defaultHash(0 / 0)).toBe(defaultHash(NaN));
    for (const k of [1, -1, 0.5, 1e300, -Infinity, 'abc', '', 10n, true, false, null, undefined, {}, Symbol()]) {
      const h = defaultHash(k);
      expect(h | 0).toBe(h);
    }
  });

  it('strings through the finalizer: the low bits differ for ids differing in one place', () => {
    const lows = new Set<number>();
    for (let i = 0; i < 32; i++) lows.add(defaultHash(`aaaaaaaaaaaaaaaaaaaaa${String.fromCharCode(65 + i)}`) & 31);
    expect(lows.size).toBeGreaterThan(16);
    expect(fmix32(0)).toBe(0);
  });

  it('integers spread over the low bits', () => {
    const lows = new Set<number>();
    for (let i = 0; i < 1024; i += 32) lows.add(defaultHash(i) & 31);
    expect(lows.size).toBeGreaterThan(16);
  });

  it('bigints by their decimal string; booleans, null and undefined apart', () => {
    expect(defaultHash(12345678901234567890n)).toBe(defaultHash(12345678901234567890n));
    const fixed = new Set([true, false, null, undefined].map(defaultHash));
    expect(fixed.size).toBe(4);
  });

  it('objects and symbols by identity; a lookup assigns no id', () => {
    const o = {};
    const f = () => {};
    const s = Symbol('x');
    const g = Symbol.for('molstar-immutable-test');
    for (const k of [o, f, s, g]) {
      expect(Number.isNaN(lookupDefaultHash(k))).toBe(true);
      const h = defaultHash(k);
      expect(defaultHash(k)).toBe(h);
      expect(lookupDefaultHash(k)).toBe(h);
    }
    expect(defaultHash({})).not.toBe(defaultHash({}));
    // a map read of an unseen object answers absent and assigns nothing
    const unseen = {};
    const m = HashMap.of(Array.from({ length: 20 }, (_, i) => [{ i }, i] as const));
    expect(m.get(unseen)).toBe(undefined);
    expect(m.has(unseen)).toBe(false);
    expect(m.delete(unseen)).toBe(m);
    expect(Number.isNaN(lookupDefaultHash(unseen))).toBe(true);
    expect(HashSet.of([o, s, g, f]).has(g)).toBe(true);
  });
});

describe('configs', () => {
  it('resolve to the defaults when empty; the same functions are the same config', () => {
    expect(resolveConfig(undefined)).toBe(DefaultConfig);
    expect(resolveConfig({})).toBe(DefaultConfig);
    expect(resolveConfig({ hash: defaultHash, equals: sameValueZero })).toBe(DefaultConfig);
    const eq = (a: unknown, b: unknown) => a === b;
    expect(sameConfig(resolveConfig({ valueEquals: eq }), resolveConfig({ valueEquals: eq }))).toBe(true);
    expect(sameConfig(resolveConfig({ valueEquals: eq }), DefaultConfig)).toBe(false);
    // valueEquals without valueHash: hashCode over keys and size
    expect(resolveConfig({ valueEquals: eq }).valueHash).toBe(undefined);
  });

  it('sameValueZero', () => {
    expect(sameValueZero(0, -0)).toBe(true);
    expect(sameValueZero(NaN, NaN)).toBe(true);
    expect(sameValueZero('a', 'a')).toBe(true);
    expect(sameValueZero({}, {})).toBe(false);
    expect(sameValueZero('NaN', NaN)).toBe(false);
  });
});
