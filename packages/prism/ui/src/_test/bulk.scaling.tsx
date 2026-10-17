// @vitest-environment jsdom
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it } from 'vitest';
import { CommandManager } from '../commands/manager.js';
import type { CommandDefinition } from '../commands/types.js';
import { flattenTree } from '../components/virtual-list/virtual-list.js';
import { render } from '../dom/render.js';
import { For, Show } from '../dom/control-flow.js';
import { flush, flushFrame, watch } from '../dom/scope.js';
import { observableState } from '../reactive/observable-state.js';
import { CommandMenu } from '../views/command-menu/command-menu.js';
import { checkLinear, quadraticControl, scaling } from './scaling.js';

// prism-ui's bulk paths grow about linearly (`test:scaling`, scaling.ts): the command menu filtering N commands, a tree
// of N siblings flattened; the quadratic control is caught.

function commands(n: number): CommandDefinition<object>[] {
  const out: CommandDefinition<object>[] = [];
  for (let i = 0; i < n; i++)
    out.push({
      id: `cmd.${i}`,
      title: `Command ${i % 7 ? 'alpha' : 'beta'} ${i}`,
      category: `C${i % 9}`,
      run: () => ({ message: '' }),
    });
  return out;
}

function manager(n: number) {
  const m = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  m.register(commands(n));
  return m;
}

it('the command manager matching a query over N commands', () => {
  const s = scaling(
    {
      setup: manager,
      run: (m) => {
        for (const q of ['al', 'alp', 'alpha 1', 'beta', 'zz']) m.query(q);
      },
    },
    1000,
  );
  checkLinear('command manager: query over N commands', s);
});

it('the command menu filtering N commands, a keystroke', () => {
  const s = scaling(
    {
      setup: (n) => {
        const box = document.createElement('div');
        document.body.append(box);
        const dispose = render(() => <CommandMenu manager={manager(n)} onClose={() => {}} />, box);
        return { box, dispose, input: box.querySelector('input')! };
      },
      run: ({ input }) => {
        for (const q of ['al', 'alp', 'alph']) {
          input.value = q;
          input.dispatchEvent(new Event('input', { bubbles: true }));
          flushFrame();
        }
      },
      teardown: ({ box, dispose }) => {
        dispose();
        box.remove();
      },
    },
    500,
  );
  checkLinear('command menu: filter N commands', s);
});

it('a tree of N siblings flattened (a reflow)', () => {
  interface N {
    readonly id: number;
    readonly children?: readonly N[];
  }
  const s = scaling(
    {
      setup: (n): N[] => [
        { id: -1, children: Array.from({ length: n }, (_, i) => ({ id: i, children: [{ id: -2 }] })) },
      ],
      run: (roots) => {
        flattenTree(
          roots,
          (x) => x.children,
          () => true,
        );
      },
    },
    5000,
  );
  checkLinear('flattenTree: N siblings', s);
});

it('N rows following one shared state, torn down together', () => {
  const s = scaling(
    {
      setup: (n) => {
        const shared = observableState(0);
        const items = observableState<readonly number[]>(Array.from({ length: n }, (_, i) => i));
        const box = document.createElement('div');
        document.body.append(box);
        const dispose = render(
          () => (
            <For each={items}>
              {(i) => {
                const el = document.createElement('div');
                watch([shared], () => (el.textContent = `${i} ${shared.value}`));
                return el;
              }}
            </For>
          ),
          box,
        );
        flush();
        return { items, box, dispose };
      },
      run: ({ items }) => {
        items.set([]);
        flush();
      },
      teardown: ({ box, dispose }) => {
        dispose();
        box.remove();
      },
    },
    2000,
  );
  checkLinear('N rows on one state, torn down', s);
});

it('one change marking N rows and their regions (watches at alternating depths)', () => {
  const s = scaling(
    {
      setup: (n) => {
        const shared = observableState(0);
        const items = observableState<readonly number[]>(Array.from({ length: n }, (_, i) => i));
        const box = document.createElement('div');
        document.body.append(box);
        const dispose = render(
          () => (
            <For each={items}>
              {(i) => {
                const el = document.createElement('div');
                const inner = document.createElement('span');
                // The row's watch, then one in a region below it: marked in turn, depth d, d + 1, d, d + 1, …
                watch([shared], () => (el.title = `${i} ${shared.value}`));
                el.append(
                  (
                    <Show when={observableState(true).asReadonly}>
                      {() => {
                        watch([shared], () => (inner.textContent = `${shared.value}`));
                        return inner;
                      }}
                    </Show>
                  ) as Node,
                );
                return el;
              }}
            </For>
          ),
          box,
        );
        flush();
        return { shared, box, dispose };
      },
      run: ({ shared }) => {
        shared.set(shared.value + 1);
        flush();
      },
      teardown: ({ box, dispose }) => {
        dispose();
        box.remove();
      },
    },
    2000,
  );
  checkLinear('one change, N watches at alternating depths', s);
});

it('a path that copies the whole list per item is caught (positive control)', () => {
  const s = scaling(quadraticControl, 2000);
  checkLinear('control: a copy of the list per item', s, false);
});
