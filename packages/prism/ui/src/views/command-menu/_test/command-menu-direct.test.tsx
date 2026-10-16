/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandMenu } from '../command-menu.js';

// The command menu opened for one command from elsewhere (`direct`): its form is the menu; one Esc closes it.

afterEach(() => {
  document.body.innerHTML = '';
});

const tick = () => new Promise((r) => setTimeout(r, 0));

function open(direct: boolean) {
  const manager = new CommandManager<object>({
    context: observableState<object>({}).asReadonly,
    formats: { get: () => undefined },
  });
  manager.register([
    {
      id: 'look.show',
      title: 'Show as',
      category: 'Look',
      params: () => ({ type: PD.Text('cartoon') }),
      run: () => ({}),
    },
  ]);
  manager.open('look.show', undefined, { type: 'surface' });
  const container = document.createElement('div');
  document.body.append(container);
  const onClose = vi.fn();
  render(() => <CommandMenu manager={manager} onClose={onClose} direct={direct} />, container);
  const input = () => container.querySelector<HTMLInputElement>('input[role=combobox]')!;
  const esc = () =>
    input().dispatchEvent(new KeyboardEvent('keydown', { key: 'Escape', bubbles: true, cancelable: true }));
  return { manager, container, onClose, esc };
}

describe('a command menu opened for one command', () => {
  it('closes with one Esc, and says so', async () => {
    const t = open(true);
    await tick();
    expect(t.container.textContent).toContain('esc close');
    t.esc();
    expect(t.manager.active.value).toBeUndefined();
    expect(t.onClose).toHaveBeenCalledTimes(1);
  });

  it('opened from the list, its Esc says it goes back to the commands', async () => {
    const t = open(false);
    await tick();
    expect(t.container.textContent).toContain('esc commands');
  });
});
