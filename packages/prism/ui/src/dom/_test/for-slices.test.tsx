/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Scheduler } from '@molstar/core/task/util/scheduler';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { observableState } from '../../reactive/observable-state.js';
import { For, timeSlicing } from '../control-flow.js';
import type { Child } from '../jsx-types.js';
import { render } from '../render.js';
import { flush } from '../scope.js';

// A long list over time slices (control-flow.ts, `For`): each slice makes only the entries the last one deferred, each
// in its place among those already placed; what ends on screen is what one reconcile gives.

let continuations: (() => void)[] = [];
let clock = 0;
beforeEach(() => {
  timeSlicing(true);
  continuations = [];
  // Each reading of the clock costs 20 ms: a slice makes its first 64 entries and is spent.
  vi.spyOn(performance, 'now').mockImplementation(() => (clock += 20));
  vi.spyOn(Scheduler, 'setImmediate').mockImplementation((fn: () => void) => {
    continuations.push(fn);
    return 0;
  });
});
afterEach(() => {
  vi.restoreAllMocks();
  timeSlicing('auto');
});

/** Runs the continuations as their tasks would, one slice at a time; how many slices ran. */
function drain(): number {
  flush();
  let ran = 0;
  while (continuations.length) {
    const fns = continuations;
    continuations = [];
    for (const f of fns) f();
    ran++;
  }
  return ran;
}

function mount(fn: () => Child) {
  const container = document.createElement('div');
  const dispose = render(fn, container);
  const shown = () => [...container.querySelectorAll('span')].map((e) => e.textContent);
  return { container, dispose, shown };
}

const range = (n: number, from = 0) => Array.from({ length: n }, (_, i) => from + i);

describe('For over time slices', () => {
  it('a long list renders in slices of new entries, in order, with their indices', () => {
    const items = observableState<readonly number[]>([]);
    const { shown, dispose } = mount(() => (
      <div>
        <For each={items}>{(i, index) => <span data-index={index.value}>{`${i}`}</span>}</For>
      </div>
    ));
    items.set(range(1000));
    flush();
    expect(shown().length).toBe(64);
    expect(drain()).toBeGreaterThan(10);
    expect(shown()).toEqual(range(1000).map(String));
    dispose();
  });

  it('new entries deferred between kept ones go in their places, the kept ones untouched', () => {
    const items = observableState<readonly number[]>(range(10));
    const { container, shown, dispose } = mount(() => (
      <div>
        <For each={items}>{(i) => <span>{`${i}`}</span>}</For>
      </div>
    ));
    drain();
    const kept = [...container.querySelectorAll('span')];
    // 100 new items before each kept one, and after the last.
    const next: number[] = [];
    for (let k = 0; k < 10; k++) next.push(...range(100, 1000 + 100 * k), k);
    next.push(...range(100, 5000));
    items.set(next);
    drain();
    expect(shown()).toEqual(next.map(String));
    const spans = [...container.querySelectorAll('span')];
    for (let k = 0; k < 10; k++) expect(spans[101 * k + 100]).toBe(kept[k]);
    dispose();
  });

  it('a list changed while some of it waits renders the latest list', () => {
    const items = observableState<readonly number[]>([]);
    const { shown, dispose } = mount(() => (
      <div>
        <For each={items}>{(i) => <span>{`${i}`}</span>}</For>
      </div>
    ));
    items.set(range(500));
    flush();
    // One slice of the first list, then another list.
    continuations.splice(0).forEach((f) => f());
    items.set(range(300, 200).reverse());
    drain();
    expect(shown()).toEqual(range(300, 200).reverse().map(String));
    items.set(range(10));
    drain();
    expect(shown()).toEqual(range(10).map(String));
    dispose();
  });

  it('a fallback shown while every entry waits goes when the first is made', () => {
    const a = observableState<readonly number[]>([]);
    const b = observableState<readonly number[]>([]);
    const { container, dispose } = mount(() => (
      <div>
        <For each={a}>{(i) => <span>{`a${i}`}</span>}</For>
        <For each={b} fallback={() => <i>none</i>}>
          {(i) => <span>{`b${i}`}</span>}
        </For>
      </div>
    ));
    expect(container.querySelector('i')).not.toBeNull();
    // The first list spends the slice: all of the second waits.
    a.set(range(100));
    b.set(range(100));
    flush();
    expect(container.querySelectorAll('span').length).toBeLessThan(200);
    drain();
    const texts = [...container.querySelectorAll('span')].map((e) => e.textContent);
    expect(texts).toEqual([...range(100).map((i) => `a${i}`), ...range(100).map((i) => `b${i}`)]);
    expect(container.querySelector('i')).toBeNull();
    dispose();
  });

  it('a later item with the key of a deferred one is not rendered', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const items = observableState<readonly { id: number }[]>([]);
    const { shown, dispose } = mount(() => (
      <div>
        <For each={items} by={(x) => x.id}>
          {(x) => <span>{`${x.id}`}</span>}
        </For>
      </div>
    ));
    const list = range(200).map((id) => ({ id }));
    items.set([...list, { id: 150 }]);
    drain();
    expect(shown()).toEqual(range(200).map(String));
    expect(warn).toHaveBeenCalled();
    dispose();
  });

  it('a list the same again after its slices does nothing', () => {
    const items = observableState<readonly number[]>([]);
    let made = 0;
    const { dispose } = mount(() => (
      <div>
        <For each={items}>
          {(i) => {
            made++;
            return <span>{`${i}`}</span>;
          }}
        </For>
      </div>
    ));
    const list = range(400);
    items.set(list);
    drain();
    expect(made).toBe(400);
    items.set([...list]);
    drain();
    expect(made).toBe(400);
    dispose();
  });
});
