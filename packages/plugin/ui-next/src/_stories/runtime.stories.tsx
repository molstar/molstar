/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { map } from 'rxjs';
import { For, Show } from '../dom/control-flow.js';
import { onCleanup } from '../dom/scope.js';
import { observableState } from '../reactive/observable-state.js';
import { story } from './story-root.js';

const meta: Meta = { title: 'ui-next/Runtime' };
export default meta;

const box = 'font: 14px system-ui, sans-serif; padding: 16px; display: grid; gap: 12px; max-width: 420px;';

/** Text, attributes and classes bound to states; nothing re-renders, bindings update the DOM. */
export const Bindings: StoryObj = {
  render: () =>
    story(() => {
      const count = observableState(0);
      const name = observableState('Mol*');
      return (
        <div style={box}>
          <label>
            Name <input value={name} onInput={(e) => name.set(e.currentTarget.value)} />
          </label>
          <p>
            Hello, <b>{name}</b>! Clicked {count} {count.asReadonly.pipe(map((c) => (c === 1 ? 'time' : 'times')))}.
          </p>
          <div style="display: flex; gap: 8px">
            <button onClick={() => count.set(count.value + 1)}>Increment</button>
            <button onClick={() => count.set(0)} disabled={count.asReadonly.pipe(map((c) => c === 0))}>
              Reset
            </button>
          </div>
        </div>
      );
    }),
};

interface Item {
  id: number;
  label: string;
}

/** A keyed list: each row keeps its DOM node and scope while its item stays; rows log their own creation. */
export const KeyedList: StoryObj = {
  render: () =>
    story(() => {
      let next = 4;
      const items = observableState<readonly Item[]>([1, 2, 3].map((id) => ({ id, label: `Item ${id}` })));
      const log = observableState<readonly string[]>([]);
      const note = (line: string) => log.set([line, ...log.value].slice(0, 6));
      const shuffle = () => items.set([...items.value].sort(() => Math.random() - 0.5));
      return (
        <div style={box}>
          <div style="display: flex; gap: 8px">
            <button onClick={() => items.set([...items.value, { id: next, label: `Item ${next++}` }])}>Add</button>
            <button onClick={shuffle}>Shuffle</button>
            <button onClick={() => items.set(items.value.slice(1))}>Remove first</button>
          </div>
          <ol data-testid="list">
            <For each={items} fallback={() => <i>No items.</i>}>
              {(item, index) => {
                note(`created ${item.label}`);
                onCleanup(() => note(`disposed ${item.label}`));
                return (
                  <li>
                    {item.label} <small>(position {index})</small>
                  </li>
                );
              }}
            </For>
          </ol>
          <pre style="margin: 0; color: #666">{log.asReadonly.pipe(map((l) => l.join('\n')))}</pre>
        </div>
      );
    }),
};

/** A branch rendered while a flag holds, disposed when it flips. */
export const Conditional: StoryObj = {
  render: () =>
    story(() => {
      const open = observableState(false);
      return (
        <div style={box}>
          <button onClick={() => open.set(!open.value)}>
            {open.asReadonly.pipe(map((o) => (o ? 'Hide details' : 'Show details')))}
          </button>
          <Show when={open} fallback={() => <i>Details are hidden.</i>}>
            {() => <p>Details: rendered in their own scope.</p>}
          </Show>
        </div>
      );
    }),
};
