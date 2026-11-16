/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  BehaviorSubject,
  batch,
  derive,
  KeyedSelector,
  KeyedState,
  RecordState,
  type ReadonlyState,
  State,
  Subject,
} from '../index.js';

// Keyed states (reactive §8), record states and named sources (reactive §4).

describe('KeyedState', () => {
  it('tells a key only when its value changes; a set of another key tells nothing', () => {
    const ks = new KeyedState<string, number>();
    const seen: unknown[] = [];
    ks.key('a').listen((v) => seen.push(['a', v]));
    ks.key('b').listen((v) => seen.push(['b', v]));
    ks.set('a', 1);
    ks.set('a', 1);
    ks.set('c', 5);
    ks.delete('a');
    ks.delete('zz');
    expect(seen).toEqual([
      ['a', 1],
      ['a', undefined],
    ]);
  });

  it('replace tells the followed keys whose value changed, came or went', () => {
    const ks = new KeyedState<string, number>();
    ks.replace(
      new Map([
        ['a', 1],
        ['b', 2],
      ]),
    );
    const seen: unknown[] = [];
    for (const k of ['a', 'b', 'c']) ks.key(k).listen((v) => seen.push([k, v]));
    ks.replace(
      new Map([
        ['a', 1],
        ['c', 3],
      ]),
    );
    expect(seen).toEqual([
      ['b', undefined],
      ['c', 3],
    ]);
  });

  it('a key is a graph node: derive follows it, and it is dropped with its last follower', () => {
    const ks = new KeyedState<string, number>((x, y) => x === y);
    const key = ks.key('a');
    const doubled = derive(key, (v) => (v ?? 0) * 2);
    expect(doubled.value).toBe(0);
    ks.set('a', 2);
    expect(doubled.value).toBe(4);
    const seen: number[] = [];
    const stop = doubled.listen((v) => seen.push(v));
    expect(ks._registry.size).toBe(1);
    ks.set('a', 3);
    stop();
    expect(ks._registry.size).toBe(0);
    expect(seen).toEqual([6]);
  });

  it('batches with the states: a key set twice in a batch is told once', () => {
    const ks = new KeyedState<string, number>();
    const seen: unknown[] = [];
    ks.key('a').listen((v) => seen.push(v));
    batch(() => {
      ks.set('a', 1);
      ks.set('a', 2);
    });
    expect(seen).toEqual([2]);
  });
});

describe('KeyedSelector', () => {
  it('a change of the current value visits only the two keys it concerns', () => {
    const current = new State('r1');
    const selector = new KeyedSelector(current);
    const flags: ReadonlyState<boolean>[] = [];
    let computes = 0;
    const told: string[] = [];
    for (let i = 0; i < 1000; i++) {
      const k = `r${i}`;
      const flag = derive(selector.key(k), (v) => {
        computes++;
        return v;
      });
      flag.listen((v) => told.push(`${k}:${v}`));
      flags.push(flag);
    }
    computes = 0;
    current.set('r500');
    expect(told).toEqual(['r1:false', 'r500:true']);
    expect(computes).toBe(2);
    expect(selector.isCurrent('r500')).toBe(true);
  });

  it('a key and another value derived from the source agree within one change', () => {
    const current = new State(1);
    const selector = new KeyedSelector(current);
    const both = derive([selector.key(2), current], (is2, c) => `${is2}:${c}`);
    const seen: string[] = [];
    both.listen((v) => seen.push(v));
    current.set(2);
    current.set(3);
    expect(seen).toEqual(['true:2', 'false:3']);
  });

  it('follows a stream without a value, only while a key is followed', () => {
    const source = new Subject<string>();
    const selector = new KeyedSelector(source);
    const a = selector.is('a');
    expect(a.value).toBe(false);
    const seen: boolean[] = [];
    const stop = a.listen((v) => seen.push(v));
    expect(source.observed).toBe(true);
    source.next('a');
    source.next('b');
    stop();
    expect(source.observed).toBe(false);
    expect(seen).toEqual([true, false]);
  });

  it('follows a BehaviorSubject', () => {
    const source = new BehaviorSubject('x');
    const selector = new KeyedSelector(source);
    const seen: boolean[] = [];
    selector.key('y').listen((v) => seen.push(v));
    source.next('y');
    expect(seen).toEqual([true]);
  });
});

describe('RecordState and named sources', () => {
  it('update sets the fields given when any differs; an equal copy tells nothing', () => {
    const r = new RecordState({ a: 1, b: 'x' });
    const seen: unknown[] = [];
    r.listen((v) => seen.push(v));
    r.update({ a: 1 });
    r.set({ a: 1, b: 'x' });
    r.update({ b: 'y' });
    expect(seen).toEqual([{ a: 1, b: 'y' }]);
  });

  it('derive over a record of sources gets a record of their values', () => {
    const a = new State(1);
    const b = new BehaviorSubject('x');
    const d = derive({ a, b }, ({ a, b }) => `${b}${a}`);
    const seen: string[] = [];
    d.listen((v) => seen.push(v));
    a.set(2);
    b.next('y');
    expect(d.value).toBe('y2');
    expect(seen).toEqual(['x2', 'y2']);
  });
});

describe('keyed states (property)', () => {
  it('every followed key ends on the map value through random sets, deletes, replaces and batches', () => {
    let seed = 3;
    const rand = (n: number) => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed % n;
    };
    const ks = new KeyedState<number, number>();
    const last = new Map<number, number | undefined>();
    const stops: (() => void)[] = [];
    for (let k = 0; k < 10; k += 2) {
      last.set(k, ks.get(k));
      stops.push(ks.key(k).listen((v) => last.set(k, v)));
    }
    const step = () => {
      const r = rand(10);
      if (r < 6) ks.set(rand(10), rand(4));
      else if (r < 8) ks.delete(rand(10));
      else ks.replace(new Map(Array.from({ length: rand(10) }, () => [rand(10), rand(4)] as [number, number])));
    };
    for (let i = 0; i < 500; i++) {
      if (rand(5) === 0) batch(() => [step(), step()]);
      else step();
      for (const [k, v] of last) expect(v).toBe(ks.get(k));
    }
    for (const stop of stops) stop();
    expect(ks._registry.size).toBe(0);
  });
});
