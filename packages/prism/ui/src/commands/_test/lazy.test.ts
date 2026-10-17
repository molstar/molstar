/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { observableState } from '../../reactive/observable-state.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition, CommandFormats } from '../types.js';

// Lazy commands (commands.md §4, `registerLazy`): made by the first read after a change, or in a microtask while the
// menu wants them; a change alone makes nothing.

const formats: CommandFormats<object> = { get: () => undefined };
const cmd = (id: string): CommandDefinition<object> => ({ id, title: id, category: 'Test', run: () => ({}) });

function setup() {
  const manager = new CommandManager<object>({ context: observableState({}).asReadonly, formats });
  let ids = ['a'];
  const make = vi.fn(() => ({ commands: ids.map(cmd), groups: [] }));
  const lazy = manager.registerLazy({ make });
  return { manager, make, lazy, setIds: (next: string[]) => (ids = next) };
}

describe('lazy commands', () => {
  it('are made by the first read, once per change, and not by a change alone', () => {
    const { manager, make, lazy, setIds } = setup();
    expect(make).not.toHaveBeenCalled();
    expect(manager.get('a')?.id).toBe('a');
    manager.query('');
    manager.children(undefined);
    expect(make).toHaveBeenCalledTimes(1);

    for (let i = 0; i < 10; i++) lazy.stale();
    setIds(['b']);
    expect(make).toHaveBeenCalledTimes(1);
    expect(manager.get('b')?.id).toBe('b');
    expect(manager.get('a')).toBeUndefined();
    expect(make).toHaveBeenCalledTimes(2);
  });

  it('keys do not make them', () => {
    const { manager, make } = setup();
    manager.handleKey({ key: 'x' });
    expect(make).not.toHaveBeenCalled();
  });

  it('while wanted, are made at once and in a microtask after each change, not by reads', async () => {
    const { manager, make, lazy, setIds } = setup();
    const release = manager.want();
    expect(make).toHaveBeenCalledTimes(1);
    expect(manager.commands.value.map((c) => c.id)).toEqual(['a']);
    setIds(['b']);
    lazy.stale();
    manager.query('');
    expect(make).toHaveBeenCalledTimes(1);
    await Promise.resolve();
    expect(make).toHaveBeenCalledTimes(2);
    expect(manager.commands.value.map((c) => c.id)).toEqual(['b']);
    release();
    lazy.stale();
    await Promise.resolve();
    expect(make).toHaveBeenCalledTimes(2);
  });

  it('a make that says nothing changed keeps what was made; dispose takes them away', () => {
    const { manager, lazy } = setup();
    const a = manager.get('a');
    const kept = manager.registerLazy({ make: () => undefined });
    lazy.stale();
    expect(manager.get('a')).not.toBe(a);
    lazy.dispose();
    expect(manager.get('a')).toBeUndefined();
    kept.dispose();
  });
});
