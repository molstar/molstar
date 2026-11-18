/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { DefaultConfig, resolveConfig, sameConfig, sameValueZero } from '../config.js';
import { defaultHash, fmix32, hashNumber, lookupDefaultHash, memoHash } from '../hash.js';
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
    const m = HashMap.of<object, number>(Array.from({ length: 20 }, (_, i) => [{ i }, i]));
    expect(m.get(unseen)).toBe(undefined);
    expect(m.has(unseen)).toBe(false);
    expect(m.delete(unseen)).toBe(m);
    expect(Number.isNaN(lookupDefaultHash(unseen))).toBe(true);
    expect(HashSet.of([o, s, g, f]).has(g)).toBe(true);
  });
});

describe('int32 keys and the string memo', () => {
  it('an int32 hashes by its value through fmix32; other numbers by their bits', () => {
    for (const n of [1, -1, 7, 2 ** 31 - 1, -(2 ** 31)]) expect(hashNumber(n)).toBe(fmix32(n));
    expect(hashNumber(-0)).toBe(0);
    expect(hashNumber(2 ** 31)).not.toBe(fmix32(2 ** 31));
    expect(hashNumber(0.5)).toBe(defaultHash(0.5));
  });

  it('memoHash gives defaultHash’s values, bounded, and assigns no id on a lookup', () => {
    const memo = memoHash(4);
    const keys = Array.from({ length: 50 }, (_, i) => `key-${i}-0123456789abcdef`);
    for (let r = 0; r < 3; r++) for (const k of keys) expect(memo(k)).toBe(defaultHash(k));
    for (const k of [1, 0.5, true, null, 10n]) expect(memo(k)).toBe(defaultHash(k));
    const a = HashMap.of(
      keys.map((k, i) => [k, i] as const),
      { hash: memo },
    );
    const b = HashMap.of(keys.map((k, i) => [k, i] as const).reverse(), { hash: memo });
    expect(a.equals(b)).toBe(true);
    expect(a.hashCode()).toBe(b.hashCode());
    for (const [k, i] of keys.map((k, i) => [k, i] as const)) expect(a.get(k)).toBe(i);
    // another memo is another config
    expect(a.equals(HashMap.of(a.toArray(), { hash: memoHash() }))).toBe(false);
    const unseen = {};
    expect(HashMap.of([[{}, 1]], { hash: memoHash() }).has(unseen)).toBe(false);
    expect(Number.isNaN(lookupDefaultHash(unseen))).toBe(true);
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
