/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { shallowEqual as shallowEqualAny } from '../index.js';
import {
  deepClone,
  deepEqual,
  isPlainRecord,
  sameData,
  sameShallow,
  shallowEqual,
  shallowEqualArrays,
  shallowEqualObjects,
} from '../object.js';

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

describe('equality', () => {
  it('isPlainRecord: an object literal or a null-prototype object only', () => {
    expect(isPlainRecord({ a: 1 })).toBe(true);
    expect(isPlainRecord(Object.create(null))).toBe(true);
    expect(isPlainRecord([])).toBe(false);
    expect(isPlainRecord(null)).toBe(false);
    expect(isPlainRecord(new (class {})())).toBe(false);
    expect(isPlainRecord(new Map())).toBe(false);
  });

  it('shallowEqualObjects: own enumerable string keys on both sides, values by Object.is', () => {
    expect(shallowEqualObjects({ a: 1, b: 'x' }, { b: 'x', a: 1 })).toBe(true);
    expect(shallowEqualObjects({ a: 1 }, { a: 1, b: undefined })).toBe(false);
    expect(shallowEqualObjects({ a: 1, b: undefined }, { a: 1 })).toBe(false);
    expect(shallowEqualObjects({ a: Number.NaN }, { a: Number.NaN })).toBe(true);
    expect(shallowEqualObjects({ a: 0 }, { a: -0 })).toBe(false);
    expect(shallowEqualObjects({ a: {} }, { a: {} })).toBe(false);
    const inherited = Object.assign(Object.create({ b: 2 }), { a: 1 });
    expect(shallowEqualObjects(inherited, { a: 1 })).toBe(true);
    const hidden = Object.defineProperty({ a: 1 }, 'b', { value: 2, enumerable: false });
    expect(shallowEqualObjects(hidden, { a: 1 })).toBe(true);
    const field = Symbol('field');
    expect(shallowEqualObjects({ a: 1, [field]: 1 }, { a: 1 })).toBe(true);
  });

  it('shallowEqualArrays and shallowEqual (object.ts: nullish values are equal)', () => {
    expect(shallowEqualArrays([1, Number.NaN], [1, Number.NaN])).toBe(true);
    expect(shallowEqualArrays([1, 2], [2, 1])).toBe(false);
    expect(shallowEqualArrays(Array(1), Array(2))).toBe(false);
    expect(shallowEqual({ a: 1 }, { a: 1 })).toBe(true);
    expect(shallowEqual([1, 2] as {}, [1, 2] as {})).toBe(true);
    expect(shallowEqual(null as unknown as {}, undefined as unknown as {})).toBe(true);
    expect(shallowEqual({ a: 1 }, undefined as unknown as {})).toBe(false);
  });

  it('shallowEqual (util/index.ts: arrays or objects; null and undefined differ)', () => {
    expect(shallowEqualAny([1, 2], [1, 2])).toBe(true);
    expect(shallowEqualAny({ a: 1 }, { a: 1 })).toBe(true);
    expect(shallowEqualAny([1], { 0: 1 })).toBe(false);
    expect(shallowEqualAny(null, undefined)).toBe(false);
    expect(shallowEqualAny(Number.NaN, Number.NaN)).toBe(true);
  });

  it('sameShallow: arrays and plain records one level deep; anything else only to itself', () => {
    const child = {};
    expect(sameShallow({ child }, { child })).toBe(true);
    expect(sameShallow({ child: {} }, { child: {} })).toBe(false);
    expect(sameShallow([1, 2], [1, 2])).toBe(true);
    expect(sameShallow([], {})).toBe(false);
    expect(sameShallow(Object.assign(Object.create(null), { a: 1 }), { a: 1 })).toBe(true);
    class Example {
      value = 1;
    }
    expect(sameShallow(new Example(), new Example())).toBe(false);
    expect(sameShallow(new Date(0), new Date(0))).toBe(false);
    expect(sameShallow(undefined, null)).toBe(false);
    expect(sameShallow(Number.NaN, Number.NaN)).toBe(true);
  });

  it('sameData: every level of plain arrays and records; anything else by identity', () => {
    expect(sameData({ a: [1, { b: 2 }] }, { a: [1, { b: 2 }] })).toBe(true);
    expect(sameData({ a: [1, { b: 2 }] }, { a: [1, { b: 3 }] })).toBe(false);
    expect(sameData({ d: new Date(0) }, { d: new Date(0) })).toBe(false);
    expect(sameData([1], [1, 2])).toBe(false);
  });

  it('deepEqual: plain data, dates and regular expressions', () => {
    expect(deepEqual({ a: [1, { b: new Date(5) }] }, { a: [1, { b: new Date(5) }] })).toBe(true);
    expect(deepEqual({ r: /x/g }, { r: /x/g })).toBe(true);
    expect(deepEqual({ a: 1 }, { a: 1, b: 2 })).toBe(false);
  });
});
