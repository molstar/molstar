/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { afterEach, describe, expect, it } from 'vitest';
import { commandItem, commandItems } from '../../../commands/items.js';
import { CommandManager } from '../../../commands/manager.js';
import { render } from '../../../dom/render.js';
import { observableState } from '../../../reactive/observable-state.js';
import { CommandButton } from '../command-button.js';

// Buttons and menu entries by command id (commands.md §7): label, key and whether it can run from the command.

const tick = () => new Promise((r) => setTimeout(r, 0));

afterEach(() => {
  document.body.innerHTML = '';
});

function setup() {
  const context = observableState({ busy: false });
  const ran: unknown[] = [];
  const manager = new CommandManager<{ busy: boolean }>({
    context: context.asReadonly,
    formats: { get: () => undefined },
  });
  const undo = manager.register([
    {
      id: 'camera.reset',
      title: 'Reset camera',
      category: 'View',
      icon: 'focus',
      shortcut: { key: 'r', alt: true },
      when: (c) => (c.busy ? 'Busy' : true),
      run: (args) => {
        ran.push(args);
        return {};
      },
    },
  ]);
  return { context, manager, ran, undo };
}

describe('CommandButton', () => {
  it('takes its icon, tooltip and key from the command, greys with the reason, runs it, and goes with it', async () => {
    const p = setup();
    const host = document.createElement('div');
    document.body.append(host);
    render(() => <CommandButton manager={p.manager} command="camera.reset" label />, host);
    const button = () => host.querySelector('button') as HTMLButtonElement;
    expect(button().textContent).toBe('Reset camera');
    expect(button().dataset.tip).toMatch(/^Reset camera \((⌥R|Alt\+R)\)$/);
    p.manager.keymap.remap('camera.reset', { key: 'h', alt: true });
    await tick();
    expect(button().dataset.tip).toMatch(/H\)$/);
    p.context.set({ busy: true });
    await tick();
    expect(button().disabled).toBe(true);
    expect(button().dataset.tip).toBe('Busy');
    p.context.set({ busy: false });
    await tick();
    button().click();
    await tick();
    expect(p.ran).toHaveLength(1);
    p.undo();
    await tick();
    expect(host.querySelector('button')).toBeNull();
  });
});

describe('commandItem', () => {
  it('makes a menu entry from the command as it is, none for a missing one', () => {
    const p = setup();
    const item = commandItem(p.manager, 'camera.reset', { label: 'Reset', hint: 'Fit everything', mac: false });
    expect(item).toMatchObject({
      id: 'camera.reset',
      label: 'Reset',
      icon: 'focus',
      shortcut: 'Alt+R',
      hint: 'Fit everything',
    });
    item!.run();
    expect(p.ran).toHaveLength(1);
    p.context.set({ busy: true });
    expect(commandItem(p.manager, 'camera.reset')?.disabled).toBe('Busy');
    expect(commandItem(p.manager, 'no.such')).toBeUndefined();
    expect(commandItems(p.manager, ['no.such', { command: 'camera.reset', id: 'reset' }]).map((i) => i.id)).toEqual([
      'reset',
    ]);
  });
});
