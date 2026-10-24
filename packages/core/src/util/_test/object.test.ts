/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { deepClone } from '../object.js';

describe('deepClone', () => {
  it('keeps typed arrays, Map and Set (perf-audit-2 01/B9)', () => {
    const src = {
      f: new Float32Array([1, 2, 3]),
      i: new Int32Array([4, 5]),
      m: new Map<string, { x: number }>([['a', { x: 1 }]]),
      s: new Set([1, 2]),
      a: [1, { b: 2 }],
    };
    const c = deepClone(src);
    expect(c.f).toBeInstanceOf(Float32Array);
    expect(Array.from(c.f)).toEqual([1, 2, 3]);
    expect(c.f).not.toBe(src.f);
    expect(c.i).toBeInstanceOf(Int32Array);
    expect(c.m).toBeInstanceOf(Map);
    expect(c.m.get('a')).toEqual({ x: 1 });
    expect(c.m.get('a')).not.toBe(src.m.get('a'));
    expect(c.s).toBeInstanceOf(Set);
    expect([...c.s]).toEqual([1, 2]);
    expect(c.a).toEqual([1, { b: 2 }]);
    expect(c.a[1]).not.toBe(src.a[1]);
    c.f[0] = 10;
    expect(src.f[0]).toBe(1);
  });

  it('copies a DataView and a typed array window', () => {
    const buf = new Uint8Array([1, 2, 3, 4]);
    const view = deepClone(new DataView(buf.buffer, 1, 2));
    expect(view.byteLength).toBe(2);
    expect(view.getUint8(0)).toBe(2);
    const win = deepClone(buf.subarray(1, 3));
    expect(Array.from(win)).toEqual([2, 3]);
    expect(win.buffer).not.toBe(buf.buffer);
  });
});
