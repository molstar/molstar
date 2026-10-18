/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { describe, expect, it, vi } from 'vitest';
import { observableState } from '../../reactive/observable-state.js';
import { chordFromEvent, chordId } from '../keys.js';
import { Keymap } from '../keymap.js';
import { CommandManager } from '../manager.js';
import type { CommandArgs, CommandDefinition, KeyBindings, KeymapStore } from '../types.js';

// The keymap's layers, keybinding scopes, groups' shared arguments and the command menu's browsing (commands.md §5–6).

type Ctx = { log: string[] };

function manager(defs: CommandDefinition<Ctx>[], options: { app?: KeyBindings; store?: KeymapStore } = {}) {
  const context = observableState<Ctx>({ log: [] });
  const preview = { set: vi.fn(), clear: vi.fn() };
  const m = new CommandManager<Ctx>({
    context: context.asReadonly,
    formats: { get: () => undefined },
    preview,
    browseDelay: 0,
    keymap: options,
  });
  m.register(defs);
  return { m, log: context.value.log, preview };
}

const cmd = (id: string, extra: Partial<CommandDefinition<Ctx>> = {}): CommandDefinition<Ctx> => ({
  id,
  title: id,
  category: 'Test',
  run: (_, ctx) => {
    ctx.log.push(id);
    return { message: id };
  },
  ...extra,
});

/** A fake element with a scope, its parent's next. */
const el = (scope?: string, parent?: object) => ({
  tagName: 'DIV',
  getAttribute: (n: string) => (n === 'data-msp-scope' ? (scope ?? null) : null),
  parentNode: parent ?? null,
});

describe('keymap layers', () => {
  it('takes the person over the app over the command, and saves the person’s', () => {
    const saved: KeyBindings[] = [];
    const store: KeymapStore = { load: () => ({ b: { key: 'y' } }), save: (b) => void saved.push(b) };
    const keymap = new Keymap({ app: { a: { key: 'x' }, c: null }, store });
    const a = { id: 'a', shortcut: { key: 'a' } };
    const b = { id: 'b', shortcut: { key: 'b' } };
    const c = { id: 'c', shortcut: { key: 'c' } };
    expect(keymap.shortcut(a)).toEqual({ key: 'x' });
    expect(keymap.shortcut(b)).toEqual({ key: 'y' });
    expect(keymap.shortcut(c)).toBeUndefined();
    keymap.remap('a', { key: 'q', mod: true });
    expect(keymap.shortcut(a)).toEqual({ key: 'q', mod: true });
    expect(keymap.isRemapped('a')).toBe(true);
    expect(keymap.defaultShortcut(a)).toEqual({ key: 'x' });
    keymap.remap('c', { key: 'z' });
    keymap.reset('a');
    expect(keymap.shortcut(a)).toEqual({ key: 'x' });
    expect(saved.at(-1)).toEqual({ b: { key: 'y' }, c: { key: 'z' } });
    keymap.resetAll();
    expect(saved.at(-1)).toEqual({});
  });

  it('drops what a store gives that is not a chord, and survives a store that throws', () => {
    const keymap = new Keymap({
      store: {
        load: () => ({ a: { key: 'k' }, b: 'nope', c: { mod: true }, d: null }) as unknown as KeyBindings,
        save: () => {
          throw new Error('full');
        },
      },
    });
    expect(keymap.layers.value.user).toEqual({ a: { key: 'k' }, d: null });
    expect(() => keymap.remap('x', { key: 'x' })).not.toThrow();
    expect(
      new Keymap({
        store: {
          load: () => {
            throw new Error('blocked');
          },
          save: () => {},
        },
      }).layers.value.user,
    ).toEqual({});
  });

  it('reports commands on one key in one scope, not a scoped one over a global one', () => {
    const keymap = new Keymap();
    const defs = [
      cmd('a', { shortcut: { key: 'k', mod: true } }),
      cmd('b', { shortcut: { key: 'K', mod: true } }),
      cmd('c', { shortcut: { key: 'k', mod: true }, scope: 'canvas' }),
      cmd('d', { shortcut: { key: '?', shift: true } }),
      cmd('e', { shortcut: { key: '?' } }),
    ];
    const conflicts = keymap.conflicts(defs);
    expect(conflicts.map((c) => [c.scope, c.commands.map((x) => x.id)])).toEqual([
      ['global', ['a', 'b']],
      ['global', ['d', 'e']],
    ]);
    expect(keymap.sharing('a', { key: 'k', mod: true }, defs).map((c) => c.id)).toEqual(['b']);
    expect(keymap.sharing('c', { key: 'k', mod: true }, defs)).toEqual([]);
  });
});

describe('chords from events', () => {
  it('reads letters, symbols, ⌥ letters by code, and refuses modifiers alone', () => {
    expect(chordFromEvent({ key: 'K', shiftKey: true, metaKey: true }, true)).toEqual({
      key: 'k',
      mod: true,
      shift: true,
    });
    expect(chordFromEvent({ key: '?', shiftKey: true }, true)).toEqual({ key: '?' });
    expect(chordFromEvent({ key: '†', code: 'KeyT', altKey: true }, true)).toEqual({ key: 't', alt: true });
    expect(chordFromEvent({ key: 'ArrowRight' }, false)).toEqual({ key: 'ArrowRight' });
    expect(chordFromEvent({ key: 'Shift', shiftKey: true }, true)).toBeUndefined();
    expect(chordFromEvent({ key: 'Tab' }, true)).toBeUndefined();
    // The platform's other modifier is never a command's.
    expect(chordFromEvent({ key: 'k', ctrlKey: true }, true)).toBeUndefined();
    expect(chordFromEvent({ key: 'k', ctrlKey: true }, false)).toEqual({ key: 'k', mod: true });
    expect(chordId({ key: '?', shift: true })).toBe(chordId({ key: '?' }));
  });
});

describe('keys through the keymap', () => {
  it('runs the remapped key, not the old one; a cleared key runs nothing', async () => {
    const { m, log } = manager(
      [cmd('a', { shortcut: { key: 'a', alt: true } }), cmd('b', { shortcut: { key: 'b', alt: true } })],
      {
        app: { b: null },
      },
    );
    expect(m.handleKey({ key: 'b', altKey: true }, false)).toBe(false);
    m.keymap.remap('a', { key: 'j', alt: true });
    expect(m.shortcut(m.get('a')!)).toEqual({ key: 'j', alt: true });
    expect(m.handleKey({ key: 'a', altKey: true }, false)).toBe(false);
    expect(m.handleKey({ key: 'j', altKey: true }, false)).toBe(true);
    await Promise.resolve();
    expect(log).toEqual(['a']);
  });
});

describe('keybinding scopes', () => {
  it('takes a scoped key only inside its scope, and there before a global one', async () => {
    const { m, log } = manager([
      cmd('global.next', { shortcut: { key: 'ArrowRight' } }),
      cmd('tree.open', { shortcut: { key: 'ArrowRight' }, scope: 'tree' }),
      cmd('canvas.only', { shortcut: { key: 'x' }, scope: 'canvas' }),
    ]);
    const panel = el('panel:structure');
    const tree = el('tree', panel);
    const leaf = el(undefined, tree);
    // One run at a time: each settles before the next key.
    const settle = () => new Promise((r) => setTimeout(r));
    expect(m.handleKey({ key: 'ArrowRight', target: leaf }, false)).toBe(true);
    await settle();
    expect(m.handleKey({ key: 'ArrowRight', target: el(undefined, panel) }, false)).toBe(true);
    await settle();
    // Not in the canvas: its key is not taken.
    expect(m.handleKey({ key: 'x', target: leaf }, false)).toBe(false);
    // A composed path (across shadow roots) is read as it is.
    expect(m.handleKey({ key: 'x', composedPath: () => [el() as never, el('canvas') as never] }, false)).toBe(true);
    await settle();
    expect(log).toEqual(['tree.open', 'global.next', 'canvas.only']);
  });

  it('lets a global command run where the scoped one on its key cannot', async () => {
    const { m, log } = manager([
      cmd('global', { shortcut: { key: 'd', mod: true } }),
      cmd('scoped', { shortcut: { key: 'd', mod: true }, scope: 'canvas', when: () => 'Nothing to do here' }),
    ]);
    expect(m.handleKey({ key: 'd', ctrlKey: true, target: el('canvas') }, false)).toBe(true);
    await Promise.resolve();
    expect(log).toEqual(['global']);
  });
});

describe('shared arguments', () => {
  const groups = [
    {
      id: 'select',
      title: 'Select',
      shared: {
        name: 'mode',
        title: 'Mode',
        options: [
          ['set', 'Replace'],
          ['add', 'Add'],
        ] as const,
        default: 'set',
      },
    },
    { id: 'select.query', title: 'By query', parent: 'select' },
  ];

  it('gives a command in the group (at any depth) the group’s value, and a form whose params take it', async () => {
    const seen: CommandArgs[] = [];
    const { m } = manager([
      cmd('select.ligand', { group: 'select.query', run: (a) => (seen.push(a), { message: '' }) }),
      cmd('select.text', {
        group: 'select',
        params: () => ({
          text: PD.Text(''),
          mode: PD.Select('add', [
            ['set', 'set'],
            ['add', 'add'],
          ]),
        }),
        run: () => ({ message: '' }),
      }),
      cmd('outside', { run: (a) => (seen.push(a), { message: '' }) }),
    ]);
    m.registerGroups(groups);
    expect(m.sharedValue('select')).toBe('set');
    await (m.open('select.ligand') as { result: Promise<unknown> }).result;
    m.setShared('select', 'add');
    await (m.open('select.ligand') as { result: Promise<unknown> }).result;
    await (m.open('outside') as { result: Promise<unknown> }).result;
    expect(seen).toEqual([{ mode: 'set' }, { mode: 'add' }, {}]);
    m.setShared('select', 'set');
    const form = m.open('select.text');
    expect(form.kind).toBe('form');
    if (form.kind === 'form') expect(form.invocation.state.value.values.mode).toBe('set');
    // A value no option has falls back to the default.
    m.setShared('select', 'nonsense');
    expect(m.sharedValue('select')).toBe('set');
  });
});

describe('browsing', () => {
  it('previews and counts the highlighted command only, and clears the preview when the highlight moves', async () => {
    const counts: AbortSignal[] = [];
    const { m, preview } = manager([
      cmd('a', {
        preview: (args) => ({ what: 'a', ...args }),
        count: (_, __, signal) => {
          counts.push(signal);
          return Promise.resolve('3 residues');
        },
      }),
      cmd('b'),
      cmd('c', { params: () => ({ x: PD.Text('') }), preview: () => 'form only' }),
    ]);
    m.browse('a');
    await Promise.resolve();
    await Promise.resolve();
    expect(preview.set).toHaveBeenCalledWith({ what: 'a' });
    expect(m.browsing.value).toEqual({ id: 'a', count: '3 residues' });
    m.browse('b');
    expect(preview.clear).toHaveBeenCalledTimes(1);
    expect(counts[0].aborted).toBe(true);
    expect(m.browsing.value).toEqual({ id: 'b' });
    // A command with params is previewed by its form, not its row.
    m.browse('c');
    expect(preview.set).toHaveBeenCalledTimes(1);
    m.browse(undefined);
    expect(m.browsing.value).toBeUndefined();
  });

  it('ends browsing when a command opens', async () => {
    const { m, preview } = manager([cmd('a', { preview: () => 'p' })]);
    m.browse('a');
    await Promise.resolve();
    expect(preview.set).toHaveBeenCalledTimes(1);
    m.open('a');
    expect(preview.clear).toHaveBeenCalled();
    expect(m.browsing.value).toBeUndefined();
  });
});
