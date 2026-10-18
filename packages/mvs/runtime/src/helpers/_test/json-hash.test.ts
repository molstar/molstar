/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { jsonHash } from '../json-hash.js';

const hash = (v: unknown) => jsonHash(v, new WeakMap());

describe('jsonHash', () => {
  it('is equal exactly when JSON.stringify is', () => {
    const values: unknown[] = [
      { a: 1, b: [1, 2] },
      { b: [1, 2], a: 1 },
      { a: 1, b: [1, 2], c: undefined },
      { a: 1, b: [1, 2], f: () => 1 },
      [undefined, 1],
      [null, 1],
      [NaN, 1],
      { x: -0 },
      { x: 0 },
      { d: new Date(0) },
      { d: '1970-01-01T00:00:00.000Z' },
      new Float32Array([1, 2]),
      { 0: 1, 1: 2 },
      [],
      {},
      '',
      'a',
      1,
      '1',
      null,
      { long: 'x'.repeat(100) },
      { long: 'x'.repeat(101) },
      { nested: { deep: { list: Array.from({ length: 50 }, (_, i) => i) } } },
      { nested: { deep: { list: Array.from({ length: 50 }, (_, i) => (i === 49 ? 0 : i)) } } },
      // large values (hashed from their parts), and their JSON-equal variants
      { big: Array.from({ length: 300 }, (_, i) => i), name: 'a' },
      { big: Array.from({ length: 300 }, (_, i) => i), name: 'a', none: undefined },
      { name: 'a', big: Array.from({ length: 300 }, (_, i) => i) },
      { big: Array.from({ length: 300 }, (_, i) => (i === 299 ? -1 : i)), name: 'a' },
      Object.fromEntries(Array.from({ length: 300 }, (_, i) => [`k${i}`, i % 3 ? 'red' : 'blue'])),
      Object.fromEntries(Array.from({ length: 300 }, (_, i) => [`k${i}`, i % 3 ? 'red' : i === 0 ? 'green' : 'blue'])),
    ];
    for (const a of values) {
      for (const b of values) {
        expect(hash(a) === hash(b), `${JSON.stringify(a)} vs ${JSON.stringify(b)}`).toBe(
          JSON.stringify(a) === JSON.stringify(b),
        );
      }
    }
  });

  it('reuses the hashes of parts it has seen', () => {
    const cache = new WeakMap<object, string>();
    const big = Array.from({ length: 1000 }, (_, i) => i);
    const h1 = jsonHash({ big, n: 1 }, cache);
    expect(cache.has(big)).toBe(true);
    expect(jsonHash({ big, n: 1 }, cache)).toBe(h1);
    expect(jsonHash({ big: big.slice(), n: 1 }, new WeakMap())).toBe(h1);
    expect(jsonHash({ big, n: 2 }, cache)).not.toBe(h1);
  });
});
