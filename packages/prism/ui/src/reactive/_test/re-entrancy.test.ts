/** @vitest-environment jsdom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Subject } from 'rxjs';
import { describe, expect, it } from 'vitest';
import { bind, createRoot, flush, requireScope, type Scope, watch } from '../../dom/scope.js';
import { keyedSelector } from '../keyed-selector.js';
import { keyedState } from '../keyed-state.js';
import { RecordState, State } from '../state.js';

// Re-entrancy (rendering.md §2.2): a set while listeners are told is told after, in order; `value` is the latest set,
// whoever reads it; every listener ends on the latest value.

describe('a set while the state tells its listeners', () => {
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

  it('a record update merges onto the latest set', () => {
    const r = new RecordState<{ a: number; b: number }>({ a: 0, b: 0 });
    let once = true;
    r.listen(() => {
      if (!once) return;
      once = false;
      r.set({ a: 1, b: 0 });
      r.update({ b: 2 });
    });
    r.set({ a: 5, b: 5 });
    expect(r.value).toEqual({ a: 1, b: 2 });
  });

  it('every listener of a keyed value ends on the latest, in order', () => {
    const k = keyedState<string, number>();
    const key = k.key('x');
    const first: (number | undefined)[] = [];
    const second: (number | undefined)[] = [];
    key.listen((v) => {
      first.push(v);
      if (v === 1) k.set('x', 2);
    });
    key.listen((v) => second.push(v));
    k.set('x', 1);
    expect(first).toEqual([1, 2]);
    expect(second).toEqual([1, 2]);
  });

  it('a key flag over a stream ends on the latest for every listener', () => {
    const picked = new Subject<string>();
    const sel = keyedSelector<string>(picked);
    const flag = sel.key('a');
    const first: boolean[] = [];
    const second: boolean[] = [];
    flag.listen((v) => {
      first.push(v);
      if (v) picked.next('b');
    });
    flag.listen((v) => second.push(v));
    picked.next('a');
    expect(first).toEqual([true, false]);
    expect(second).toEqual([true, false]);
  });
});

describe('a key flag over a stream that replays', () => {
  it('the replayed current value is not a change: a watch renders once on start', () => {
    const picked = new BehaviorSubject('a');
    const flag = keyedSelector<string>(picked).key('a');
    let renders = 0;
    const dispose = createRoot((d) => {
      watch([flag], () => {
        renders++;
        void flag.value;
      });
      return d;
    });
    flush();
    expect(renders).toBe(1);
    picked.next('b');
    flush();
    expect(renders).toBe(2);
    dispose();
  });
});

describe('bind', () => {
  it('a set made during its first write is written', () => {
    const s = new State(0);
    const written: number[] = [];
    const dispose = createRoot((d) => {
      bind(s, (v) => {
        if (v === 0) s.set(1);
        written.push(v);
      });
      return d;
    });
    expect(written[written.length - 1]).toBe(1);
    expect(s.value).toBe(1);
    dispose();
  });

  it('a set made during the write on resume is written', () => {
    const s = new State(0);
    const written: number[] = [];
    let scope!: Scope;
    const dispose = createRoot((d) => {
      scope = requireScope('test');
      bind(s, (v) => {
        if (v === 5) s.set(6);
        written.push(v);
      });
      return d;
    });
    scope.setActive(false);
    s.set(5);
    scope.setActive(true);
    expect(written[written.length - 1]).toBe(6);
    dispose();
  });
});
