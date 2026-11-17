/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { ReactiveConfig } from '../config.js';
import {
  BehaviorSubject,
  batch,
  derive,
  distinctUntilChanged,
  isState,
  listen,
  map,
  type ReadonlyState,
  State,
  Subject,
  toState,
} from '../index.js';

afterEach(() => {
  ReactiveConfig.onUnhandledError = null;
});

// The rules of states (reactive §4): ordered delivery, batching, re-entrancy; and what push-dirty, pull-values gives:
// a diamond told once, a derived value nobody reads not computed, a chain that stops at an unchanged value.

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

  it('value is the latest set: a listener counting from it counts every step', () => {
    const c = new State(0);
    let once = true;
    c.listen(() => {
      if (!once) return;
      once = false;
      c.set(c.value + 1);
      c.set(c.value + 1);
    });
    c.set(10);
    expect(c.value).toBe(12);
  });

  it('skips equal values (Object.is, or its equals)', () => {
    const s = new State(Number.NaN);
    const seen: number[] = [];
    s.listen((v) => seen.push(v));
    s.set(Number.NaN);
    const r = new State({ a: 1 }, (x, y) => x.a === y.a);
    const told: unknown[] = [];
    r.listen((v) => told.push(v));
    r.set({ a: 1 });
    r.set({ a: 2 });
    expect(seen).toEqual([]);
    expect(told).toEqual([{ a: 2 }]);
  });

  it('is an Observable: the current value, then each change; a set during the replay is told after it', () => {
    const s = new State(1);
    const seen: string[] = [];
    s.subscribe((v) => {
      seen.push(`a${v}`);
      if (v === 1) s.set(2);
    });
    s.set(3);
    expect(seen).toEqual(['a1', 'a2', 'a3']);
    const piped: number[] = [];
    s.pipe(map((v) => v * 10)).subscribe((v) => piped.push(v));
    expect(piped).toEqual([30]);
  });

  it('asReadonly: a view with value, listen and subscribe, no set', () => {
    const s = new State(1);
    const view = s.asReadonly;
    expect(s.asReadonly).toBe(view);
    expect('set' in view).toBe(false);
    const seen: number[] = [];
    view.listen((v) => seen.push(v));
    s.set(2);
    expect(view.value).toBe(2);
    expect(seen).toEqual([2]);
    expect(isState(view)).toBe(true);
  });

  it('reports what a listener throws, and the others are told', () => {
    const reported: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => reported.push(e);
    const s = new State(0);
    const seen: number[] = [];
    s.listen(() => {
      throw new Error('bad');
    });
    s.listen((v) => seen.push(v));
    s.set(1);
    expect(reported).toEqual([new Error('bad')]);
    expect(seen).toEqual([1]);
  });
});

describe('batch', () => {
  it('tells each state once with its final value, in the order first set; reads inside see the new values', () => {
    const a = new State(0);
    const b = new State(0);
    const seen: string[] = [];
    a.listen((v) => seen.push(`a${v}`));
    b.listen((v) => seen.push(`b${v}`));
    const result = batch(() => {
      a.set(1);
      b.set(1);
      a.set(2);
      expect(a.value).toBe(2);
      return 'r';
    });
    expect(result).toBe('r');
    expect(seen).toEqual(['a2', 'b1']);
  });

  it('tells nothing for a state set back to where it started, and holds until the outermost batch ends', () => {
    const a = new State(0);
    const seen: number[] = [];
    a.listen((v) => seen.push(v));
    batch(() => {
      a.set(1);
      a.set(0);
    });
    batch(() => {
      batch(() => a.set(5));
      expect(seen).toEqual([]);
    });
    expect(seen).toEqual([5]);
  });

  it('tells after a throwing batch', () => {
    const a = new State(0);
    const seen: number[] = [];
    a.listen((v) => seen.push(v));
    expect(() =>
      batch(() => {
        a.set(1);
        throw new Error('x');
      }),
    ).toThrow();
    expect(seen).toEqual([1]);
  });
});

describe('derive', () => {
  it('is lazy: computes on read, follows nothing while nothing listens, caches until something changes', () => {
    const a = new State(1);
    let runs = 0;
    const d = derive(a, (x) => {
      runs++;
      return x * 2;
    });
    expect(runs).toBe(0);
    expect(d.value).toBe(2);
    expect(d.value).toBe(2);
    expect(runs).toBe(1);
    expect(a.observed).toBe(false);
    a.set(2);
    expect(d.value).toBe(4);
    expect(runs).toBe(2);
  });

  it('is hot and shared while listened: one computation per change, equal results not told', () => {
    const a = new State(1);
    let runs = 0;
    const parity = derive(a, (x) => {
      runs++;
      return x % 2;
    });
    const seen: string[] = [];
    parity.listen((v) => seen.push(`x${v}`));
    parity.listen((v) => seen.push(`y${v}`));
    a.set(3);
    a.set(4);
    expect(seen).toEqual(['x0', 'y0']);
    expect(runs).toBe(3);
    // Followed by the derived value (a graph dependent), not by a listener of its own.
    expect(a.observed).toBe(true);
  });

  it('stops following its sources when the last listener goes', () => {
    const a = new State(1);
    const d = derive(a, (x) => x + 1);
    const e = derive(d, (x) => x + 1);
    const stop = e.listen(() => {});
    expect(a.observed).toBe(true);
    expect((a as any)._refs).toBe(1);
    stop();
    expect(a.observed || d.observed).toBe(false);
    a.set(5);
    expect(e.value).toBe(7);
  });

  it('derives from several states, with a spread past four', () => {
    const s = [1, 2, 3, 4, 5].map((v) => new State(v));
    const sum = derive(s, (...xs) => xs.reduce((a, b) => a + b, 0));
    const pair = derive([s[0], s[1]], (x, y) => `${x}-${y}`);
    const seen: number[] = [];
    sum.listen((v) => seen.push(v));
    s[4].set(10);
    expect(seen).toEqual([20]);
    expect(pair.value).toBe('1-2');
  });

  it('a diamond computes and tells once, never with a stale side', () => {
    const a = new State(1);
    const double = derive(a, (x) => x * 2);
    const triple = derive(a, (x) => x * 3);
    let runs = 0;
    const sum = derive([double, triple, a], (x, y, z) => {
      runs++;
      return `${z}:${x + y}`;
    });
    const seen: string[] = [];
    sum.listen((v) => seen.push(v));
    runs = 0;
    a.set(2);
    expect(seen).toEqual(['2:10']);
    expect(runs).toBe(1);
  });

  it('a chain stops at the first unchanged value', () => {
    const a = new State(1);
    const sign = derive(a, (x) => Math.sign(x));
    let downstream = 0;
    const label = derive(sign, (s) => {
      downstream++;
      return s > 0 ? 'positive' : 'not';
    });
    const seen: string[] = [];
    label.listen((v) => seen.push(v));
    downstream = 0;
    a.set(5);
    a.set(7);
    expect(downstream).toBe(0);
    a.set(-1);
    expect(seen).toEqual(['not']);
    expect(downstream).toBe(1);
  });

  it('a derived value nothing reads is not computed when its source changes', () => {
    const a = new State(1);
    let runs = 0;
    const unread = derive(a, (x) => {
      runs++;
      return x;
    });
    const listened = derive(a, (x) => x);
    listened.listen(() => {});
    a.set(2);
    a.set(3);
    expect(runs).toBe(0);
    expect(unread.value).toBe(3);
  });

  it('follows a BehaviorSubject directly (no toState)', () => {
    const b = new BehaviorSubject(1);
    const d = derive([b, new State(10)], (x, y) => x + y);
    const seen: number[] = [];
    const stop = d.listen((v) => seen.push(v));
    b.next(2);
    b.next(2);
    expect(seen).toEqual([12]);
    stop();
    expect(b.observed).toBe(false);
    b.next(5);
    expect(d.value).toBe(15);
  });

  it('a throwing compute throws on read, is reported in a delivery, and recovers', () => {
    const reported: unknown[] = [];
    ReactiveConfig.onUnhandledError = (e) => reported.push(e);
    const a = new State(1);
    const d = derive(a, (x) => {
      if (x === 2) throw new Error('two');
      return x;
    });
    const seen: number[] = [];
    d.listen((v) => seen.push(v));
    a.set(2);
    expect(() => d.value).toThrow('two');
    a.set(3);
    expect(reported).toEqual([new Error('two')]);
    expect(seen).toEqual([3]);
  });

  it('a subscriber of a derived value gets its value, then each change', () => {
    const a = new State(1);
    const d = derive(a, (x) => x * 2);
    const seen: number[] = [];
    d.pipe(distinctUntilChanged()).subscribe((v) => seen.push(v));
    a.set(2);
    expect(seen).toEqual([2, 4]);
  });
});

describe('toState and listen', () => {
  it('toState feeds a state from a stream, once, until disposed', () => {
    const s = new Subject<number>();
    const st = toState(s, 0);
    const seen: number[] = [];
    st.listen((v) => seen.push(v));
    s.next(1);
    s.next(1);
    st.dispose();
    s.next(2);
    expect(seen).toEqual([1]);
    expect(s.observed).toBe(false);
  });

  it('listen: a replayed current value is no change', () => {
    const b = new BehaviorSubject(1);
    let n = 0;
    const stop = listen(b.asObservable(), () => n++);
    b.next(2);
    stop();
    b.next(3);
    expect(n).toBe(1);
    const st = new State(0);
    const stopSt = listen(st, () => n++);
    st.set(1);
    stopSt();
    expect(n).toBe(2);
  });
});

describe('random graphs (property)', () => {
  it("every listener ends on a from-scratch recompute of its node's value", () => {
    let seed = 11;
    const rand = (n: number) => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed % n;
    };
    for (let round = 0; round < 40; round++) {
      const states = Array.from({ length: 4 }, (_, i) => new State(i));
      type Node = { state: ReadonlyState<number>; compute: () => number };
      const nodes: Node[] = states.map((s) => ({ state: s, compute: () => s.value }));
      for (let i = 0; i < 12; i++) {
        const picks = Array.from({ length: 1 + rand(3) }, () => nodes[rand(nodes.length)]);
        const k = rand(5);
        const fn = (...xs: number[]) => xs.reduce((a, b) => a + b, k) % 7;
        const state =
          picks.length === 1
            ? derive(picks[0].state, (x) => fn(x))
            : derive(
                picks.map((p) => p.state),
                fn,
              );
        nodes.push({ state, compute: () => fn(...picks.map((p) => p.compute())) });
      }
      const last = new Map<Node, number>();
      const stops: (() => void)[] = [];
      for (const node of nodes) {
        if (rand(2) === 0) continue;
        last.set(node, node.state.value);
        stops.push(node.state.listen((v) => last.set(node, v)));
      }
      for (let step = 0; step < 20; step++) {
        const set = () => states[rand(states.length)].set(rand(10));
        if (rand(4) === 0) batch(() => [set(), set(), set()]);
        else set();
        for (const [node, value] of last) expect(value).toBe(node.compute());
        for (const node of nodes) expect(node.state.value).toBe(node.compute());
      }
      for (const stop of stops) stop();
      for (const s of states) expect(s.observed).toBe(false);
    }
  });
});

describe('set back in a batch', () => {
  it('restores the version: a derived value does not recompute, nothing is told', () => {
    const a = new State(1);
    let runs = 0;
    const d = derive(a, (x) => {
      runs++;
      return x * 2;
    });
    const seen: number[] = [];
    d.listen((v) => seen.push(v));
    a.listen((v) => seen.push(v));
    runs = 0;
    batch(() => {
      a.set(2);
      a.set(1);
    });
    expect(runs).toBe(0);
    expect(seen).toEqual([]);
    a.set(3);
    expect(seen).toEqual([3, 6]);
  });

  it('a new value after a set-back gets a version of its own (a cold read in between is not stale)', () => {
    const a = new State(1);
    const d = derive(a, (x) => x * 10);
    batch(() => {
      a.set(2);
      expect(d.value).toBe(20);
      a.set(1);
      a.set(3);
      expect(d.value).toBe(30);
    });
  });
});
