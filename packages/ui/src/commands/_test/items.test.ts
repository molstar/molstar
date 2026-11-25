/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { commandItem, commandItems } from '../items.js';
import { chordKeys, formatChord, KeyPlatform } from '../keys.js';
import { CommandManager } from '../manager.js';

// Menu entries by command id and chord display (commands §6.5, §7).

afterEach(() => {
  KeyPlatform.mac = false;
});

function setup() {
  const context = new State({ busy: false, log: [] as string[] });
  const manager = new CommandManager<typeof context.value>({ context, formats: { get: () => undefined } });
  manager.register([
    {
      id: 'view.reset',
      title: 'Reset camera',
      category: 'View',
      icon: 'refresh-cw',
      shortcut: { key: 'r', mod: true },
      when: (c) => (c.busy ? 'Busy.' : true),
      run: (_a, c) => (c.log.push('reset'), {}),
    },
    {
      id: 'look.show',
      title: 'Show…',
      category: 'Look',
      params: () => ({ what: PD.Text('') }),
      run: () => ({}),
    },
  ]);
  return { manager, context };
}

describe('commandItem', () => {
  it('takes the title, icon, key after the keymap and the reason it cannot run; runs through the manager', async () => {
    const { manager, context } = setup();
    const item = commandItem(manager, 'view.reset', { mac: true });
    expect(item).toMatchObject({ id: 'view.reset', label: 'Reset camera', icon: 'refresh-cw', shortcut: '⌘R' });
    expect(item?.disabled).toBeUndefined();
    item?.run();
    await Promise.resolve();
    expect(context.value.log).toEqual(['reset']);
    manager.keymap.remap('view.reset', { key: 'x', alt: true });
    context.set({ ...context.value, busy: true });
    const blocked = vi.fn();
    const again = commandItem(manager, 'view.reset', { label: 'Reset', mac: false, onBlocked: blocked });
    expect(again).toMatchObject({ label: 'Reset', shortcut: 'Alt+X', disabled: 'Busy.' });
    again?.run();
    expect(blocked).toHaveBeenCalledWith('Busy.');
  });

  it('opens the invocation of a command with params at the values given; none for an unknown id', () => {
    const { manager } = setup();
    commandItem(manager, 'look.show', { values: { what: 'ligand' } })?.run();
    expect(manager.active.value?.values()).toEqual({ what: 'ligand' });
    expect(commandItem(manager, 'nope')).toBeUndefined();
    expect(commandItems(manager, ['nope', 'view.reset', { command: 'look.show', id: 'x' }]).map((i) => i.id)).toEqual([
      'view.reset',
      'x',
    ]);
  });
});

describe('chord display', () => {
  it('follows the platform the host set', () => {
    expect(formatChord({ key: 'k', mod: true })).toBe('Ctrl+K');
    KeyPlatform.mac = true;
    expect(formatChord({ key: 'k', mod: true })).toBe('⌘K');
    expect(chordKeys({ key: 'Backspace', alt: true, shift: true })).toEqual(['⌥', '⇧', '⌫']);
  });
});
