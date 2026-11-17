/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  BindingPresets,
  bindingTable,
  DefaultActions,
  GenomeBrowserBindings,
  type ModifierState,
  modifiersOf,
  SelectionBindings,
} from '../bindings.js';
import { CommandChangeKinds, changeKindOf, isMethodAction } from '../commands.js';

const none = {};
const shift = { shift: true };
const mod = { mod: true };
const alt = { alt: true };
const space = { space: true };

describe('gesture bindings', () => {
  it('the default preset resolves each gesture by the hit’s target', () => {
    const t = bindingTable();
    expect(t.issues).toEqual([]);
    expect(t.resolve('hover', none, 'cell')).toBe('hover');
    expect(t.resolve('hover', none, 'empty')).toBe('hover');
    expect(t.resolve('drag', none, 'cell')).toBeUndefined();
    expect(t.resolve('drag', none, 'axis')).toBe('pan');
    expect(t.resolve('drag', space, 'cell')).toBe('pan');
    expect(t.resolve('drag', alt, 'axis')).toBeUndefined();
    expect(t.resolve('wheel', none, 'item')).toBe('pan');
    expect(t.resolve('wheel', shift, 'item')).toBe('pan');
    expect(t.resolve('wheel', mod, 'cell')).toBe('zoom');
    expect(t.resolve('pinch', none, 'cell')).toBe('zoom');
    expect(t.resolve('double-click', none, 'item')).toBe('view.focus');
    expect(t.resolve('double-click', none, 'empty')).toBe('view.fit');
    expect(t.resolve('click', none, 'cell')).toBeUndefined();
  });

  it('modifiers do not stop the hover: it is bound with any modifiers, tried after the exact ones', () => {
    const t = bindingTable();
    const held: ModifierState[] = [none, shift, mod, alt, space, { shift: true, mod: true, alt: true, space: true }];
    for (const m of held)
      for (const target of ['cell', 'item', 'empty'] as const) expect(t.resolve('hover', m, target)).toBe('hover');
    // an exact binding beats the any-modifiers one
    const o = bindingTable(undefined, [[{ gesture: 'hover', modifiers: ['alt'], target: 'any', action: 'pan' }]]);
    expect(o.resolve('hover', alt, 'cell')).toBe('pan');
    expect(o.resolve('hover', shift, 'cell')).toBe('hover');
  });

  it('a preset where a drag on data pans, as genome browsers do', () => {
    const t = bindingTable(GenomeBrowserBindings);
    expect(t.issues).toEqual([]);
    expect(t.resolve('drag', none, 'cell')).toBe('pan');
    expect(Object.keys(BindingPresets)).toEqual(['default', 'genome-browser']);
  });

  it('overrides replace a binding of the same key, null removes one; conflicts and missing actions are reported', () => {
    const t = bindingTable(undefined, [
      [
        { gesture: 'double-click', target: 'item', action: 'view.fit' },
        { gesture: 'wheel', modifiers: ['mod'], target: 'any', action: null },
      ],
    ]);
    expect(t.resolve('double-click', none, 'item')).toBe('view.fit');
    expect(t.resolve('double-click', none, 'cell')).toBe('view.focus');
    expect(t.resolve('wheel', mod, 'cell')).toBeUndefined();

    const bad = bindingTable(
      [
        { gesture: 'click', target: 'cell', action: 'view.fit' },
        { gesture: 'click', target: 'cell', action: 'view.focus' },
        { gesture: 'click', target: 'item', action: 'no.such' },
        { gesture: 'key', target: 'any', action: 'view.fit' },
      ],
      [],
    );
    expect(bad.issues.map((i) => i.code)).toEqual([
      'seq.viewer.binding.conflict',
      'seq.viewer.binding.action',
      'seq.viewer.binding.key',
    ]);
    // an action that does not exist is left out: the gesture does nothing
    expect(bad.resolve('click', none, 'item')).toBeUndefined();
    // the selection's bindings wait for its model and commands
    expect(bindingTable(SelectionBindings).issues.every((i) => i.code === 'seq.viewer.binding.action')).toBe(true);
    // an application's own command exists when it says so, with its kind
    const own = bindingTable([{ gesture: 'click', target: 'item', action: 'app.open' }], [], {
      actions: [...DefaultActions, 'app.open'],
      kinds: { 'app.open': 'session' },
    });
    expect(own.issues).toEqual([]);
    expect(own.changeKind('click', none, 'item')).toBe('session');
  });

  it('modifiers are platform-aware: Cmd on macOS, Ctrl elsewhere; one record reused', () => {
    expect(modifiersOf({ metaKey: true }, true).mod).toBe(true);
    expect(modifiersOf({ ctrlKey: true }, true).mod).toBe(false);
    expect(modifiersOf({ ctrlKey: true }, false).mod).toBe(true);
    const out: ModifierState = {};
    expect(modifiersOf({ shiftKey: true }, false, true, out)).toBe(out);
    expect(out).toEqual({ shift: true, mod: false, alt: false, space: true });
  });
});

describe('commands’ change kinds', () => {
  it('every action of the presets is a method or a command with a kind', () => {
    for (const preset of Object.values(BindingPresets))
      for (const b of preset)
        expect(isMethodAction(b.action!) || changeKindOf(b.action!) !== undefined, b.action!).toBe(true);
    const t = bindingTable();
    expect(t.changeKind('double-click', none, 'empty')).toBe('navigate');
    expect(t.changeKind('wheel', none, 'cell')).toBeUndefined();
    expect(changeKindOf('nothing')).toBeUndefined();
    expect(changeKindOf('view.fit', { 'view.fit': 'session' })).toBe('session');
    expect(Object.values(CommandChangeKinds).every((k) => ['patch', 'session', 'navigate'].includes(k))).toBe(true);
  });
});
