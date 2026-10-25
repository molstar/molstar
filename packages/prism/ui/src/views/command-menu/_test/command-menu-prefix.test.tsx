/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it, vi } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import type { CommandMenuPrefix } from '../../../commands/types.js';
import { render } from '../../../dom/render.js';
import { flushFrame } from '../../../dom/scope.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu } from '../command-menu.js';

// A prefix's rows (commands.md §5): a word prefix with its space, a row that runs, and a completion that puts its text
// in the menu and keeps it open.

afterEach(() => {
  document.body.innerHTML = '';
});

const tick = () => new Promise((r) => setTimeout(r, 0));

function open(prefix: CommandMenuPrefix) {
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  const container = document.createElement('div');
  document.body.append(container);
  const onClose = vi.fn();
  render(() => <CommandMenu manager={manager} onClose={onClose} prefixes={[prefix]} />, container);
  const input = () => container.querySelector<HTMLInputElement>('input[role=combobox]')!;
  const type = (text: string) => {
    input().value = text;
    input().dispatchEvent(new Event('input', { bubbles: true }));
    // The list is filtered once a frame (rendering.md §2.3).
    flushFrame();
  };
  const key = (k: string) =>
    input().dispatchEvent(new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true }));
  const rows = () => [...container.querySelectorAll('[role=option]')].map((r) => r.textContent);
  return { container, onClose, input, type, key, rows };
}

describe("a command-menu prefix's rows", () => {
  it('fills the text from a completion and keeps the menu open; a run closes it', async () => {
    const ran = vi.fn();
    const t = open({
      char: 'say ',
      title: 'Say',
      items: (text) =>
        text.startsWith('hello')
          ? [{ id: 'run', label: `Run: ${text}`, run: () => void ran(text) }]
          : [{ id: 'fill', label: 'hello', run: () => ({ text: 'say hello world' }) }],
    });
    await tick();
    t.type('say he');
    await tick();
    expect(t.rows()).toEqual(['hello']);
    t.key('Enter');
    await tick();
    expect(t.input().value).toBe('say hello world');
    expect(t.onClose).not.toHaveBeenCalled();
    expect(t.rows()).toEqual(['Run: hello world']);
    t.key('Enter');
    await tick();
    expect(ran).toHaveBeenCalledWith('hello world');
    expect(t.onClose).toHaveBeenCalledTimes(1);
  });
});
