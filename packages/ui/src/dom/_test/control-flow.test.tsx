/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { type ReadonlyState, State } from '@molstar/core/reactive';
import {
  activity,
  Dynamic,
  For,
  flush,
  Match,
  Motion,
  onCleanup,
  render,
  Show,
  Switch,
  timeSlicing,
  ViewScope,
} from '../index.js';

// Control flow (dom §5).

afterEach(() => {
  Motion.close = undefined;
  Motion.swap = undefined;
  timeSlicing('auto');
});

function mount(fn: () => Node, parent?: ViewScope) {
  const host = document.createElement('div');
  const dispose = render(fn, host, parent);
  return { host, dispose };
}

describe('Show', () => {
  it('renders the content while truthy, else the fallback, anew only when the value changes identity', () => {
    const when = new State<{ name: string } | undefined>(undefined);
    let renders = 0;
    const { host, dispose } = mount(() => (
      <Show when={when} fallback={() => 'none'}>
        {(v) => {
          renders++;
          return <b>{v.name}</b>;
        }}
      </Show>
    ));
    expect(host.textContent).toBe('none');
    const value = { name: 'a' };
    when.set(value);
    expect(host.textContent).toBe('a');
    when.set(value);
    when.set({ name: 'b' });
    expect(host.textContent).toBe('b');
    expect(renders).toBe(2);
    when.set(undefined);
    expect(host.textContent).toBe('none');
    dispose();
  });

  it('keeps its content across a pause; disposes the content scope when it goes', () => {
    const when = new State(true);
    const parent = new ViewScope();
    let cleanups = 0;
    let renders = 0;
    const { host } = mount(
      () => (
        <Show when={when}>
          {() => {
            renders++;
            onCleanup(() => cleanups++);
            return <i>on</i>;
          }}
        </Show>
      ),
      parent,
    );
    parent.setActive(false);
    parent.setActive(true);
    expect(renders).toBe(1);
    when.set(false);
    expect(cleanups).toBe(1);
    expect(host.textContent).toBe('');
  });

  it('reveal marks opening content; a closing plays through Motion.close before it goes', async () => {
    let finish!: () => void;
    const finished = new Promise<void>((r) => (finish = r));
    Motion.close = vi.fn(() => ({ finished, cancel() {} }) as unknown as Animation);
    const when = new State(false);
    const { host } = mount(() => (
      <Show when={when} reveal="side">
        {() => <p>x</p>}
      </Show>
    ));
    when.set(true);
    expect(host.querySelector('p')!.getAttribute('data-reveal')).toBe('side');
    when.set(false);
    expect(Motion.close).toHaveBeenCalledWith(host.querySelector('p'), 'side');
    expect(host.querySelector('p')!.hasAttribute('inert')).toBe(true);
    finish();
    await finished;
    await Promise.resolve();
    expect(host.querySelector('p')).toBeNull();
  });
});

describe('Switch, Dynamic', () => {
  it('Switch renders the first truthy case, else the fallback', () => {
    const a = new State(false);
    const b = new State<string | undefined>('b');
    const { host } = mount(() => (
      <Switch fallback={() => 'none'}>
        <Match when={a}>{() => 'A'}</Match>
        <Match when={b}>{(v) => `B:${v}`}</Match>
      </Switch>
    ));
    expect(host.textContent).toBe('B:b');
    a.set(true);
    expect(host.textContent).toBe('A');
    a.set(false);
    b.set(undefined);
    expect(host.textContent).toBe('none');
  });

  it('Dynamic renders the current component with the props', () => {
    const One = (p: { label: string }) => <b>1{p.label}</b>;
    const Two = (p: { label: string }) => <i>2{p.label}</i>;
    const c = new State<((p: { label: string }) => Node) | undefined>(One);
    const { host } = mount(() => <Dynamic component={c} label="x" />);
    expect(host.innerHTML).toContain('<b>1x</b>');
    c.set(Two);
    expect(host.innerHTML).toContain('<i>2x</i>');
    c.set(undefined);
    expect(host.textContent).toBe('');
  });
});

describe('For', () => {
  it('keeps an entry while its key and item stay; moves, adds and removes the rest', () => {
    const a = { id: 'a' };
    const b = { id: 'b' };
    const c = { id: 'c' };
    const items = new State([a, b, c]);
    let made = 0;
    const { host } = mount(() => (
      <ul>
        <For each={items} by={(x) => x.id} fallback={() => <li>empty</li>}>
          {(x) => {
            made++;
            return <li>{x.id}</li>;
          }}
        </For>
      </ul>
    ));
    const li = () => [...host.querySelectorAll('li')];
    const [la, lb, lc] = li();
    items.set([c, a, b]);
    expect(li().map((e) => e.textContent)).toEqual(['c', 'a', 'b']);
    expect(li()).toEqual([lc, la, lb]);
    expect(made).toBe(3);
    items.set([c, { id: 'a' }]);
    expect(li().map((e) => e.textContent)).toEqual(['c', 'a']);
    expect(li()[1]).not.toBe(la);
    expect(made).toBe(4);
    items.set([]);
    expect(li().map((e) => e.textContent)).toEqual(['empty']);
    items.set([b]);
    expect(li().map((e) => e.textContent)).toEqual(['b']);
  });

  it('gives each entry its index as a state; warns once for a duplicate key', () => {
    const warn = vi.spyOn(console, 'warn').mockImplementation(() => {});
    const items = new State(['x', 'y']);
    const indices: ReadonlyState<number>[] = [];
    mount(() => (
      <For each={items}>
        {(x, i) => {
          indices.push(i);
          return <span>{x}</span>;
        }}
      </For>
    ));
    items.set(['y', 'x', 'x', 'x']);
    expect(indices.map((i) => i.value)).toEqual([1, 0]);
    expect(warn).toHaveBeenCalledTimes(1);
    warn.mockRestore();
  });

  it('time slicing: a long list renders over slices, each making only the entries deferred', async () => {
    timeSlicing(true);
    const items = new State(Array.from({ length: 2000 }, (_, i) => i));
    const { host } = mount(() => (
      <div>
        <For each={items}>{(x) => <p>{x}</p>}</For>
      </div>
    ));
    const before = host.querySelectorAll('p').length;
    expect(before).toBeGreaterThan(0);
    for (let i = 0; i < 200 && host.querySelectorAll('p').length < 2000; i++) {
      await new Promise((r) => setTimeout(r, 0));
    }
    const all = [...host.querySelectorAll('p')].map((p) => Number(p.textContent));
    expect(all.length).toBe(2000);
    expect(all.every((v, i) => v === i)).toBe(true);
  });
});

describe('activity', () => {
  it('pauses its content and makes its top elements inert while inactive; the DOM stays', () => {
    const active = new State(true);
    const label = new State('a');
    const { host } = mount(() => (
      <div>
        {activity(active, () => (
          <section>{label}</section>
        ))}
      </div>
    ));
    const section = host.querySelector('section')!;
    active.set(false);
    expect(section.hasAttribute('inert')).toBe(true);
    label.set('b');
    flush();
    expect(section.textContent).toBe('a');
    expect(label.observed).toBe(false);
    active.set(true);
    expect(section.hasAttribute('inert')).toBe(false);
    expect(section.textContent).toBe('b');
    expect(host.querySelector('section')).toBe(section);
  });
});
