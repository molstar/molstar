/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '../../reactive/observable-state.js';
import { Dynamic, For, Match, Show, Switch, timeSlicing } from '../control-flow.js';
import type { Child } from '../jsx-types.js';
import { render } from '../render.js';
import { onCleanup } from '../scope.js';

function mount(fn: () => Child) {
  const container = document.createElement('div');
  const dispose = render(fn, container);
  return { container, dispose, html: () => container.innerHTML.replace(/<!--[^>]*-->/g, '') };
}

describe('Show', () => {
  it('renders the branch for a truthy value and the fallback otherwise, disposing what it leaves', () => {
    const user = observableState<{ name: string } | undefined>(undefined);
    const cleaned: string[] = [];
    const { html } = mount(() => (
      <Show when={user} fallback={() => <i>nobody</i>}>
        {(u) => {
          onCleanup(() => cleaned.push(u.name));
          return <b>{u.name}</b>;
        }}
      </Show>
    ));
    expect(html()).toBe('<i>nobody</i>');
    user.set({ name: 'a' });
    expect(html()).toBe('<b>a</b>');
    user.set({ name: 'b' });
    expect(html()).toBe('<b>b</b>');
    expect(cleaned).toEqual(['a']);
    user.set(undefined);
    expect(html()).toBe('<i>nobody</i>');
    expect(cleaned).toEqual(['a', 'b']);
  });

  it('keeps the branch while a flag stays truthy', () => {
    const on = observableState(true);
    let renders = 0;
    mount(() => (
      <Show when={on}>
        {() => {
          renders++;
          return <b />;
        }}
      </Show>
    ));
    on.set(true);
    expect(renders).toBe(1);
  });
});

describe('For', () => {
  interface Row {
    id: number;
    label: string;
  }
  const row = (id: number, label = `r${id}`): Row => ({ id, label });

  it('keeps each item’s nodes and scope while its key stays; moves, adds and removes the others', () => {
    const a = row(1);
    const b = row(2);
    const c = row(3);
    const items = observableState<readonly Row[]>([a, b, c]);
    const created: number[] = [];
    const disposed: number[] = [];
    const { container, html } = mount(() => (
      <ul>
        <For each={items}>
          {(r, index) => {
            created.push(r.id);
            onCleanup(() => disposed.push(r.id));
            return (
              <li>
                {r.label}:{index}
              </li>
            );
          }}
        </For>
      </ul>
    ));
    expect(html()).toBe('<ul><li>r1:0</li><li>r2:1</li><li>r3:2</li></ul>');
    const [liA, liB] = Array.from(container.querySelectorAll('li'));
    const d = row(4);
    items.set([c, a, d]);
    expect(html()).toBe('<ul><li>r3:0</li><li>r1:1</li><li>r4:2</li></ul>');
    expect(container.querySelectorAll('li')[1]).toBe(liA);
    expect(liB.isConnected).toBe(false);
    expect(created).toEqual([1, 2, 3, 4]);
    expect(disposed).toEqual([2]);
  });

  it('renders an item again when its key stays but the value changes', () => {
    const items = observableState<readonly Row[]>([row(1, 'old')]);
    const { container, html } = mount(() => (
      <ul>
        <For each={items} key={(r) => r.id}>
          {(r) => <li>{r.label}</li>}
        </For>
      </ul>
    ));
    const before = container.querySelector('li');
    items.set([row(1, 'new')]);
    expect(html()).toBe('<ul><li>new</li></ul>');
    expect(container.querySelector('li')).not.toBe(before);
  });

  it('keeps items with several nodes together and shows a fallback when empty', () => {
    const items = observableState<readonly string[]>(['a', 'b']);
    const { html } = mount(() => (
      <div>
        <For each={items} fallback={() => <i>empty</i>}>
          {(x) => (
            <>
              <dt>{x}</dt>
              <dd>{x.toUpperCase()}</dd>
            </>
          )}
        </For>
      </div>
    ));
    expect(html()).toBe('<div><dt>a</dt><dd>A</dd><dt>b</dt><dd>B</dd></div>');
    items.set(['b', 'a']);
    expect(html()).toBe('<div><dt>b</dt><dd>B</dd><dt>a</dt><dd>A</dd></div>');
    items.set([]);
    expect(html()).toBe('<div><i>empty</i></div>');
    items.set(['c']);
    expect(html()).toBe('<div><dt>c</dt><dd>C</dd></div>');
  });
});

describe('For, reordering', () => {
  it('ends in the list order after shuffles, inserts and removals, keeping the nodes of entries that stay', () => {
    // A fixed pseudo-random sequence.
    let seed = 7;
    const rand = (n: number) => {
      seed = (seed * 16807) % 2147483647;
      return seed % n;
    };
    const items = observableState<number[]>(Array.from({ length: 40 }, (_, i) => i));
    const { container, dispose } = mount(() => (
      <ul>
        <For each={items}>{(i) => <li>{i}</li>}</For>
      </ul>
    ));
    let next = 40;
    for (let round = 0; round < 200; round++) {
      const before = new Map([...container.querySelectorAll('li')].map((li) => [Number(li.textContent), li]));
      const list = items.value.slice();
      const op = rand(4);
      if (op === 0) {
        // Swap two.
        const a = rand(list.length);
        const b = rand(list.length);
        [list[a], list[b]] = [list[b], list[a]];
      } else if (op === 1) {
        // Shuffle.
        for (let i = list.length - 1; i > 0; i--) {
          const j = rand(i + 1);
          [list[i], list[j]] = [list[j], list[i]];
        }
      } else if (op === 2) {
        list.splice(rand(list.length + 1), 0, next++, next++);
      } else if (list.length > 2) {
        list.splice(rand(list.length), 2);
      }
      items.set(list);
      const after = [...container.querySelectorAll('li')];
      expect(after.map((li) => Number(li.textContent))).toEqual(list);
      for (const li of after) {
        const kept = before.get(Number(li.textContent));
        if (kept) expect(li).toBe(kept);
      }
    }
    dispose();
  });
});

describe('For, a long list', () => {
  it('renders over several slices once the slice is spent, in order, and follows the latest list', async () => {
    // jsdom renders nothing, so slicing is off there unless asked for.
    timeSlicing(true);
    // Each reading of the clock costs 1 ms: the 5 ms slice is spent after a few new entries past the first 64.
    let clock = 0;
    const now = jest.spyOn(performance, 'now').mockImplementation(() => (clock += 1));
    const items = observableState(Array.from({ length: 300 }, (_, i) => i));
    const { container, dispose } = mount(() => (
      <div>
        <For each={items}>{(i) => <span>{i}</span>}</For>
      </div>
    ));
    try {
      const shown = () => [...container.querySelectorAll('span')].map((e) => Number(e.textContent));
      const first = shown().length;
      expect(first).toBeGreaterThanOrEqual(64);
      expect(first).toBeLessThan(300);
      // The list changes while some of it waits: the rest follows the latest list.
      items.set(items.value.filter((i) => i % 2 === 0));
      // The continuations are fresh tasks (Mol*'s setImmediate): let them run.
      for (let n = 0; n < 500 && shown().length < 150; n++) await new Promise((r) => setTimeout(r, 0));
      expect(shown()).toEqual(items.value);
    } finally {
      dispose();
      now.mockRestore();
      timeSlicing('auto');
    }
  });
});

describe('Dynamic and Switch', () => {
  it('renders the component that is current, with the props given', () => {
    const A = (p: { n: number }) => <a>{p.n}</a>;
    const B = (p: { n: number }) => <b>{p.n}</b>;
    const current = observableState<typeof A | undefined>(A);
    const { html } = mount(() => <Dynamic component={current} n={1} />);
    expect(html()).toBe('<a>1</a>');
    current.set(B);
    expect(html()).toBe('<b>1</b>');
    current.set(undefined);
    expect(html()).toBe('');
  });

  it('renders the first case that holds', () => {
    const error = observableState<string | undefined>(undefined);
    const loading = observableState(true);
    const second = mount(() => (
      <Switch fallback={() => <i>ready</i>}>
        <Match when={loading}>{() => <i>loading</i>}</Match>
        <Match when={error}>{(m) => <i>{m}</i>}</Match>
      </Switch>
    ));
    expect(second.html()).toBe('<i>loading</i>');
    loading.set(false);
    error.set('failed');
    expect(second.html()).toBe('<i>failed</i>');
    error.set(undefined);
    expect(second.html()).toBe('<i>ready</i>');
  });
});
