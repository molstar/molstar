/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { utf8Encode } from '@molstar/core/util/utf8';
import { describe, expect, it } from 'vitest';
import { parseXyz } from '../xyz/parser.js';
import { plainResult } from './plain-result.js';

const XyzString = `3
water, ÅÉ 水
O   0.000000  0.000000  0.117300
H\t0.000000 0.757200 -0.469200
  H   0.000000  -0.757200  -0.469200
2
second
C 1.5 -2.25e1 3
N   4  5   6
`;

describe('xyz reader', () => {
  it('basic', async () => {
    const parsed = await parseXyz(XyzString).run();
    if (parsed.isError) throw new Error(parsed.message);
    const { molecules } = parsed.result;
    expect(molecules.length).toBe(2);
    const [a, b] = molecules;
    expect(a.count).toBe(3);
    expect(a.comment).toBe('water, ÅÉ 水');
    expect(a.type_symbol.toArray()).toEqual(['O', 'H', 'H']);
    expect(Array.from(a.y.toArray())).toEqual([0, 0.7572, -0.7572]);
    expect(Array.from(a.z.toArray())).toEqual([0.1173, -0.4692, -0.4692]);
    expect(b.comment).toBe('second');
    expect(b.type_symbol.toArray()).toEqual(['C', 'N']);
    expect(Array.from(b.x.toArray())).toEqual([1.5, 4]);
    expect(Array.from(b.y.toArray())).toEqual([-22.5, 5]);
  });

  it('parses the same from a string and from its bytes', async () => {
    for (const text of [XyzString, XyzString.replace(/\n/g, '\r\n')]) {
      const fromString = await parseXyz(text).run();
      const fromBytes = await parseXyz(utf8Encode(text)).run();
      expect(plainResult(fromBytes)).toEqual(plainResult(fromString));
    }
  });

  it('missing coordinates are NaN', async () => {
    const parsed = await parseXyz('1\n\nHe 1.0\n').run();
    if (parsed.isError) throw new Error(parsed.message);
    const [m] = parsed.result.molecules;
    expect(m.type_symbol.value(0)).toBe('He');
    expect(m.x.value(0)).toBe(1);
    expect(m.y.value(0)).toBeNaN();
  });
});
