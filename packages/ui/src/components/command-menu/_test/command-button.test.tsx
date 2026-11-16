/** @vitest-environment jsdom */
/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { CommandManager } from '../../../commands/manager.js';
import { KeyPlatform } from '../../../commands/keys.js';
import { flush, render } from '../../../dom/index.js';
import { CommandButton } from '../../button/command-button.js';

// A button by command id (commands §7): its tooltip with the key after the keymap, disabled with the reason, gone when
// unregistered.

const mac = KeyPlatform.mac;
afterEach(() => {
  KeyPlatform.mac = mac;
  document.body.innerHTML = '';
});

describe('CommandButton', () => {
  it('follows the command: its key, its `when`, its registration', () => {
    KeyPlatform.mac = false;
    const ctx = new State({ loaded: false });
    const run = vi.fn(() => ({}));
    const m = new CommandManager<{ loaded: boolean }>({ context: ctx.asReadonly, formats: { get: () => undefined } });
    const unregister = m.register([
      {
        id: 'view.reset',
        title: 'Reset…',
        category: 'View',
        icon: 'rotate-ccw',
        shortcut: { key: 'r', mod: true },
        when: (c) => (c.loaded ? true : 'Nothing is loaded.'),
        run,
      },
    ]);
    const container = document.createElement('div');
    document.body.append(container);
    render(() => <CommandButton manager={m} command="view.reset" label />, container);
    const button = () => container.querySelector('button');
    expect(button()?.textContent).toBe('Reset…');
    expect(button()?.getAttribute('data-tip')).toBe('Nothing is loaded.');
    expect(button()?.disabled).toBe(true);
    ctx.set({ loaded: true });
    flush();
    expect(button()?.getAttribute('data-tip')).toBe('Reset (Ctrl+R)');
    expect(button()?.disabled).toBe(false);
    m.keymap.remap('view.reset', { key: 'x', alt: true });
    flush();
    expect(button()?.getAttribute('data-tip')).toBe('Reset (Alt+X)');
    button()?.click();
    expect(run).toHaveBeenCalledTimes(1);
    unregister();
    flush();
    expect(button()).toBeNull();
  });
});
