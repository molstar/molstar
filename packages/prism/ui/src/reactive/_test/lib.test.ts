/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { createDisposableGroup } from '../disposable-group.js';
import { sameShallow } from '../same-shallow.js';
import { shallowEqual } from '../shallow-equal.js';

describe('sameShallow', () => {
  it('compares primitives by Object.is and nested values by identity', () => {
    expect(sameShallow(Number.NaN, Number.NaN)).toBe(true);
    expect(sameShallow(0, -0)).toBe(false);
    expect(sameShallow(undefined, null)).toBe(false);
    const child = {};
    expect(sameShallow({ child }, { child })).toBe(true);
    expect(sameShallow({ child: {} }, { child: {} })).toBe(false);
    expect(sameShallow(Object.assign(Object.create(null), { a: 1 }), { a: 1 })).toBe(true);
  });

  it('compares arrays, sparse lengths, and enumerable symbol fields', () => {
    expect(sameShallow([1, 2], [1, 2])).toBe(true);
    expect(sameShallow([1, 2], [2, 1])).toBe(false);
    expect(sameShallow(Array(1), Array(2))).toBe(false);
    expect(sameShallow([], {})).toBe(false);
    const field = Symbol('state');
    expect(sameShallow({ [field]: 1 }, { [field]: 1 })).toBe(true);
    expect(sameShallow({ [field]: 1 }, { [field]: 2 })).toBe(false);
  });

  it('does not collapse class instances, dates, or maps with matching enumerable fields', () => {
    class Example {
      value = 1;
    }
    expect(sameShallow(new Example(), new Example())).toBe(false);
    expect(sameShallow(new Date(0), new Date(0))).toBe(false);
    expect(sameShallow(new Map(), new Map())).toBe(false);
    const instance = new Example();
    expect(sameShallow(instance, instance)).toBe(true);
  });
});

describe('shallowEqual', () => {
  it('compares keys with Object.is', () => {
    const left = { a: 1, b: 'two' };
    const right = { a: 1, b: 'two' };
    expect(shallowEqual(left, right)).toBe(true);
    expect(shallowEqual(left, { a: 1, b: 'three' })).toBe(false);
  });

  it('counts own enumerable keys only, on both sides', () => {
    expect(shallowEqual({ a: 1 }, { a: 1, b: undefined } as { a: number })).toBe(false);
    expect(shallowEqual({ a: 1, b: undefined } as { a: number }, { a: 1 })).toBe(false);
    expect(shallowEqual({ a: Number.NaN }, { a: Number.NaN })).toBe(true);
    expect(shallowEqual({ a: 0 }, { a: -0 })).toBe(false);
    const inherited = Object.assign(Object.create({ b: 2 }), { a: 1 });
    expect(shallowEqual(inherited, { a: 1 })).toBe(true);
    expect(shallowEqual({ a: 1, b: 2 }, inherited)).toBe(false);
    const hidden = Object.defineProperty({ a: 1 }, 'b', { value: 2, enumerable: false });
    expect(shallowEqual(hidden, { a: 1 })).toBe(true);
    const field = Symbol('state');
    expect(shallowEqual({ a: 1, [field]: 1 }, { a: 1 })).toBe(false);
    expect(shallowEqual({ a: 1 }, { a: 1, [field]: 1 })).toBe(false);
    expect(shallowEqual({ a: 1, [field]: 1 }, { a: 1, [field]: 1 })).toBe(true);
  });
});

describe('createDisposableGroup', () => {
  it('disposes immediately when added to an already disposed group', () => {
    const group = createDisposableGroup();
    group.dispose();
    let disposed = false;
    group.registerDispose(() => {
      disposed = true;
    });
    expect(disposed).toBe(true);
  });

  it('removes a dynamic child without disturbing retained siblings', () => {
    const group = createDisposableGroup();
    const order: string[] = [];
    group.registerDispose(() => {
      order.push('retained');
    });
    const removeChild = group.registerDispose(() => {
      order.push('child');
    });

    removeChild();
    group.dispose();
    removeChild();

    expect(order).toEqual(['child', 'retained']);
  });

  it('attempts every cleanup in reverse order and reports the first failure', () => {
    const group = createDisposableGroup();
    const order: string[] = [];
    group.registerDispose(() => {
      order.push('first');
      throw new Error('first failure');
    });
    group.registerDispose(() => {
      order.push('second');
      throw new Error('second failure');
    });
    group.registerDispose(() => {
      order.push('third');
    });

    expect(() => group.dispose()).toThrow('second failure');
    expect(order).toEqual(['third', 'second', 'first']);
    expect(() => group.dispose()).not.toThrow();
  });
});
