/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// prism-ui's hot paths in jsdom, shared by their quick cases (hot.quick.ts) and heap guards (hot.heap.ts): a `For`
// given the same rows in a new array, typing in the command menu. jsdom's own DOM allocates: the baseline is a
// recorded constant, not zero; the heap guards hold that nothing is kept.

import { CommandManager } from '../../commands/manager.js';
import type { CommandDefinition } from '../../commands/types.js';
import { observableState } from '../../reactive/observable-state.js';
import { CommandMenu } from '../../views/command-menu/command-menu.js';
import { For } from '../control-flow.js';
import { render } from '../render.js';
import { flush, flushFrame } from '../scope.js';

export interface Case {
  step: (i: number) => void;
  sink: () => number;
}

/** A `For` of 200 rows given the same rows in another array (made ahead): nothing to do. */
export function forUnchanged(): Case {
  const rows = Array.from({ length: 200 }, (_, i) => ({ id: i, label: `Row ${i}` }));
  const a = rows.slice();
  const b = rows.slice();
  const list = observableState<readonly { id: number; label: string }[]>(a, { compare: Object.is });
  const box = document.createElement('div');
  render(
    () => (
      <For each={list} by={(r) => r.id}>
        {(r) => <div>{r.label}</div>}
      </For>
    ),
    box,
  );
  return {
    step: (i) => {
      list.set(i % 2 ? a : b);
      flush();
    },
    sink: () => box.childElementCount,
  };
}

/** Typing in the command menu over 300 commands: a letter typed and erased, the same rows with new marks. */
export function commandMenuTyping(): Case {
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  const commands: CommandDefinition<object>[] = [];
  for (let i = 0; i < 300; i++)
    commands.push({
      id: `cmd.${i}`,
      title: `Command ${i % 7 ? 'alpha' : 'beta'} ${i}`,
      category: `C${i % 9}`,
      run: () => ({ message: '' }),
    });
  manager.register(commands);
  const box = document.createElement('div');
  document.body.append(box);
  render(() => <CommandMenu manager={manager} onClose={() => {}} />, box);
  const input = box.querySelector('input')!;
  // A query narrowed by a letter and back: the same commands match, their marks move (typing).
  const texts = ['alp', 'alph'];
  return {
    step: (i) => {
      input.value = texts[i & 1];
      input.dispatchEvent(new Event('input', { bubbles: true }));
      // The frame the input asks for (the list is filtered once a frame).
      flushFrame();
    },
    sink: () => box.querySelectorAll('[role=option]').length,
  };
}
