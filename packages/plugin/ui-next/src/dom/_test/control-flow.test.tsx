/** @jest-environment jsdom */
/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { observableState } from '../../reactive/observable-state.js';
import { Dynamic, For, Match, Show, Switch } from '../control-flow.js';
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
