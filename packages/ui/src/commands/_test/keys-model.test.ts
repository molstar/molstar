/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { State } from '@molstar/core/reactive';
import { afterEach, describe, expect, it } from 'vitest';
import { KeyPlatform } from '../keys.js';
import { KeysModel } from '../keys-model.js';
import { CommandManager } from '../manager.js';
import type { CommandDefinition, KeyBindings } from '../types.js';

// Changing keys headless (commands §6.6): what Help lists, conflicts, a key recorded, cleared, reset; what each said.

const mac = KeyPlatform.mac;
afterEach(() => {
  KeyPlatform.mac = mac;
});

function setup() {
  KeyPlatform.mac = true;
  const saved: KeyBindings[] = [];
  const manager = new CommandManager<object>({
    context: new State<object>({}).asReadonly,
    formats: { get: () => undefined },
    keymap: { app: { 'app.cleared': null }, store: { load: () => undefined, save: (b) => void saved.push(b) } },
  });
  const cmd = (
    id: string,
    title: string,
    extra: Partial<CommandDefinition<object>> = {},
  ): CommandDefinition<object> => ({
    id,
    title,
    category: 'Test',
    run: () => ({}),
    ...extra,
  });
  manager.register([
    cmd('next', 'Next step', { shortcut: { key: 'ArrowRight' } }),
    cmd('reset', 'Reset camera', { shortcut: { key: 'r', alt: true } }),
    cmd('app.cleared', 'Cleared by the app', { shortcut: { key: 'x', alt: true } }),
    cmd('plain', 'No key at first'),
  ]);
  const keys = new KeysModel(manager);
  const titles = () => keys.listed.value.map((c) => c.title);
  return { manager, keys, saved, titles };
}

describe('the keys model', () => {
  it('lists the effective keys; the app’s removal is gone, and what has none can get one', () => {
    const { keys, titles } = setup();
    expect(titles()).toEqual(['Next step', 'Reset camera']);
    expect(keys.keyless.value.map((c) => c.id)).toEqual(['app.cleared', 'plain']);
    expect(keys.conflicts.value).toEqual([]);
  });

  it('records a key, says a conflict, resets and clears', () => {
    const { keys, manager, saved, titles } = setup();
    keys.setEditing(true);
    keys.record('reset');
    expect(keys.recording.value).toBe('reset');
    keys.recorded({ key: 'ArrowRight' });
    expect(keys.recording.value).toBeUndefined();
    expect(saved.at(-1)).toEqual({ reset: { key: 'ArrowRight' } });
    expect(keys.said.value).toBe('Reset camera: →, also used by Next step.');
    expect(keys.conflicts.value.map((c) => c.commands.map((x) => x.id))).toEqual([['next', 'reset']]);
    keys.reset('reset');
    expect(manager.shortcut(manager.get('reset')!)).toEqual({ key: 'r', alt: true });
    expect(keys.said.value).toBe('Reset camera: ⌥R again.');
    expect(keys.conflicts.value).toEqual([]);
    keys.clear('next');
    expect(manager.shortcut(manager.get('next')!)).toBeUndefined();
    // Cleared by the person: still listed, to reset it.
    expect(titles()).toContain('Next step');
    expect(keys.said.value).toBe('Next step: no key.');
  });

  it('gives a key to a command without one, and cancels a recording', () => {
    const { keys, manager, titles } = setup();
    keys.record('plain');
    // Listed while recorded, to show the recorder.
    expect(titles()).toContain('No key at first');
    keys.cancel();
    expect(titles()).not.toContain('No key at first');
    keys.record('plain');
    keys.recorded({ key: 'p', mod: true });
    expect(manager.shortcut(manager.get('plain')!)).toEqual({ key: 'p', mod: true });
    keys.resetAll();
    expect(manager.shortcut(manager.get('plain')!)).toBeUndefined();
    expect(keys.said.value).toBe('Every key is back to its default.');
    keys.record('nothing');
    expect(keys.recording.value).toBeUndefined();
  });

  it('says the same change again: the text differs by a trailing no-break space, so it is told', () => {
    const { keys } = setup();
    const told: string[] = [];
    const stop = keys.said.listen!((t) => told.push(t));
    keys.resetAll();
    keys.resetAll();
    keys.resetAll();
    stop();
    expect(told).toEqual([
      'Every key is back to its default.',
      'Every key is back to its default.\u00a0',
      'Every key is back to its default.',
    ]);
  });

  it('says it in the app’s words', () => {
    const { manager } = setup();
    const keys = new KeysModel(manager, {
      labels: { cleared: (c) => `${c} : aucune touche.`, resetAll: 'Tout est rétabli.' },
    });
    keys.clear('next');
    expect(keys.said.value).toBe('Next step : aucune touche.');
    keys.resetAll();
    expect(keys.said.value).toBe('Tout est rétabli.');
  });
});
