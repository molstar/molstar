/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { CommandManager, type CommandManagerOptions } from '../../../commands/manager.js';
import type { CommandDefinition } from '../../../commands/types.js';
import { flush, flushFrame, render } from '../../../dom/index.js';
import { CommandMenu, type CommandMenuProps } from '../command-menu.js';

export const tick = () => new Promise((r) => setTimeout(r, 0));

/** A manager over an empty context with no formats unless given. */
export function manager(
  commands: readonly CommandDefinition<object>[] = [],
  options: Partial<CommandManagerOptions<object>> = {},
): CommandManager<object> {
  const m = new CommandManager<object>({
    context: new State<object>({}).asReadonly,
    formats: { get: () => undefined },
    ...options,
  });
  m.register(commands);
  return m;
}

/** The menu rendered in the document, with helpers that act as the person does. */
export function mount(props: Omit<CommandMenuProps, 'onClose'> & { onClose?: () => void }) {
  const container = document.createElement('div');
  document.body.append(container);
  let closed = 0;
  let dispose: (() => void) | undefined;
  const onClose = () => {
    closed++;
    props.onClose?.();
  };
  dispose = render(() => <CommandMenu {...props} onClose={onClose} />, container);
  const input = () => container.querySelector<HTMLInputElement>('input[role=combobox]')!;
  const key = (k: string, init: KeyboardEventInit = {}, target: Element = input()) => {
    const e = new KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true, ...init });
    target.dispatchEvent(e);
    flush();
    return e;
  };
  const type = (text: string) => {
    input().value = text;
    input().dispatchEvent(new Event('input', { bubbles: true }));
    flushFrame();
    flush();
  };
  const options = () => [...container.querySelectorAll<HTMLElement>('[role=option]')];
  const rows = () => options().map((r) => r.textContent ?? '');
  return {
    container,
    input,
    key,
    type,
    options,
    rows,
    closed: () => closed,
    dispose: () => {
      dispose?.();
      dispose = undefined;
    },
  };
}
