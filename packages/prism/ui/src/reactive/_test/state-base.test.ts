/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { map, Subject } from 'rxjs';
import { describe, expect, it, vi } from 'vitest';
import { batch, derive, State, toState } from '../state.js';

// The reactive base (rendering.md §2.2): a state is its latest value plus listeners; `derive` is lazy and hot; rxjs
// only at the boundary (a state is an Observable when subscribed; `toState` feeds a state from a stream).

describe('State', () => {
  it('tells listeners each value in order, a set inside a delivery after it', () => {
    const s = new State(0);
    const seen: string[] = [];
    s.listen((v) => {
      seen.push(`a${v}`);
      if (v === 1) s.set(2);
    });
    s.listen((v) => seen.push(`b${v}`));
    s.set(1);
    expect(seen).toEqual(['a1', 'b1', 'a2', 'b2']);
    expect(s.value).toBe(2);
  });

  it('a listener removed while others are told is not told again; one added is told the next value', () => {
    const s = new State(0);
    const seen: string[] = [];
    let stopB: () => void = () => {};
    s.listen((v) => {
      seen.push(`a${v}`);
      stopB();
      if (v === 1) s.listen((w) => seen.push(`c${w}`));
    });
    stopB = s.listen((v) => seen.push(`b${v}`));
    s.set(1);
    s.set(2);
    expect(seen).toEqual(['a1', 'a2', 'c2']);
  });

  it('holds sets in a batch and tells the final value once', () => {
    const s = new State(0);
    const fn = vi.fn();
    s.listen(fn);
    batch(() => {
      s.set(1);
      s.set(2);
      expect(s.value).toBe(2);
    });
    expect(fn.mock.calls).toEqual([[2]]);
  });

  it('is an Observable at the boundary: the current value, then each new one', () => {
    const s = new State(1);
    const seen: number[] = [];
    const sub = s.asReadonly.pipe(map((v) => v * 10)).subscribe((v) => seen.push(v));
    s.set(2);
    sub.unsubscribe();
    s.set(3);
    expect(seen).toEqual([10, 20]);
  });
});

describe('derive', () => {
  it('is lazy: computes on read and follows nothing while nothing listens', () => {
    const a = new State(1);
    const fn = vi.fn((v: number) => v * 2);
    const d = derive(a, fn);
    expect(d.value).toBe(2);
    a.set(2);
    expect(d.value).toBe(4);
    const stop = d.listen!(() => {});
    a.set(3);
    expect(d.value).toBe(6);
    stop();
    const calls = fn.mock.calls.length;
    a.set(4);
    expect(fn.mock.calls.length).toBe(calls);
  });

  it('is hot and shared while listened: one computation per change, equal values skipped', () => {
    const a = new State(1);
    const fn = vi.fn((v: number) => v % 2);
    const d = derive(a, fn);
    const one = vi.fn();
    const two = vi.fn();
    d.listen!(one);
    d.listen!(two);
    fn.mockClear();
    a.set(3);
    expect(fn).toHaveBeenCalledTimes(1);
    expect(one).not.toHaveBeenCalled();
    a.set(4);
    expect(one.mock.calls).toEqual([[0]]);
    expect(two.mock.calls).toEqual([[0]]);
  });

  it('derives from several states', () => {
    const a = new State(1);
    const b = new State('x');
    const d = derive([a, b], (n, s) => `${s}${n}`);
    const seen: string[] = [];
    d.listen!((v) => seen.push(v));
    a.set(2);
    b.set('y');
    expect(seen).toEqual(['x2', 'y2']);
    expect(d.value).toBe('y2');
  });
});

describe('toState', () => {
  it('feeds a state from a stream, once, until disposed', () => {
    const events = new Subject<number>();
    const s = toState(events, 0);
    const seen: number[] = [];
    s.listen!((v) => seen.push(v));
    events.next(1);
    events.next(1);
    events.next(2);
    expect(seen).toEqual([1, 2]);
    expect(events.observed).toBe(true);
    s.dispose();
    expect(events.observed).toBe(false);
  });
});
