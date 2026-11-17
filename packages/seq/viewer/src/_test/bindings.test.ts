/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { BindingPresets, bindingTable, GenomeBrowserBindings, modifiersOf } from '../bindings.js';
import { changeKindOf, CommandChangeKinds, isMethodAction } from '../commands.js';

const none = {};
const shift = { shift: true };
const mod = { mod: true };
const space = { space: true };

describe('gesture bindings', () => {
  it('the default preset resolves each gesture by the hit’s target', () => {
    const t = bindingTable();
    expect(t.issues).toEqual([]);
    expect(t.resolve('hover', none, 'cell')).toBe('hover');
    expect(t.resolve('hover', none, 'empty')).toBe('hover');
    expect(t.resolve('click', none, 'cell')).toBe('selection.set');
    expect(t.resolve('click', none, 'item')).toBe('selection.set');
    expect(t.resolve('click', none, 'row')).toBe('selection.set');
    expect(t.resolve('click', none, 'empty')).toBe('selection.clear');
    expect(t.resolve('click', shift, 'cell')).toBe('selection.extend');
    expect(t.resolve('click', mod, 'item')).toBe('selection.toggle');
    expect(t.resolve('click', { shift: true, mod: true }, 'item')).toBeUndefined();
    expect(t.resolve('drag', none, 'cell')).toBe('drag-select');
    expect(t.resolve('drag', none, 'axis')).toBe('pan');
    expect(t.resolve('drag', space, 'cell')).toBe('pan');
    expect(t.resolve('wheel', none, 'item')).toBe('pan');
    expect(t.resolve('wheel', mod, 'cell')).toBe('zoom');
    expect(t.resolve('pinch', none, 'cell')).toBe('zoom');
    expect(t.resolve('double-click', none, 'item')).toBe('view.focus');
    expect(t.resolve('double-click', none, 'empty')).toBe('view.fit');
  });

  it('a preset where a drag on data pans, as genome browsers do', () => {
    const t = bindingTable(GenomeBrowserBindings);
    expect(t.issues).toEqual([]);
    expect(t.resolve('drag', none, 'cell')).toBe('pan');
    expect(t.resolve('drag', shift, 'item')).toBe('drag-select');
    expect(Object.keys(BindingPresets)).toEqual(['default', 'genome-browser']);
  });

  it('overrides replace a binding of the same key, null removes one; conflicts and unknown actions are reported', () => {
    const t = bindingTable(undefined, [
      [
        { gesture: 'double-click', target: 'item', action: 'selection.set' },
        { gesture: 'wheel', modifiers: ['mod'], target: 'any', action: null },
      ],
    ]);
    expect(t.resolve('double-click', none, 'item')).toBe('selection.set');
    expect(t.resolve('double-click', none, 'cell')).toBe('view.focus');
    expect(t.resolve('wheel', mod, 'cell')).toBeUndefined();

    const bad = bindingTable(
      [
        { gesture: 'click', target: 'cell', action: 'selection.set' },
        { gesture: 'click', target: 'cell', action: 'view.fit' },
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
    // an application's own command is known when its kind is given
    const own = bindingTable([{ gesture: 'click', target: 'item', action: 'app.open' }], [], {
      commands: { 'app.open': 'session' },
    });
    expect(own.issues).toEqual([]);
    expect(own.changeKind('click', none, 'item')).toBe('session');
  });

  it('modifiers are platform-aware: Cmd on macOS, Ctrl elsewhere', () => {
    expect(modifiersOf({ metaKey: true }, true).mod).toBe(true);
    expect(modifiersOf({ ctrlKey: true }, true).mod).toBe(false);
    expect(modifiersOf({ ctrlKey: true }, false).mod).toBe(true);
    expect(modifiersOf({ shiftKey: true }, false, true)).toEqual({ shift: true, mod: false, alt: false, space: true });
  });
});

describe('commands’ change kinds', () => {
  it('every action of the presets is a method or a command with a kind', () => {
    for (const preset of Object.values(BindingPresets))
      for (const b of preset)
        expect(isMethodAction(b.action!) || changeKindOf(b.action!) !== undefined, b.action!).toBe(true);
    const t = bindingTable();
    expect(t.changeKind('click', none, 'cell')).toBe('session');
    expect(t.changeKind('double-click', none, 'empty')).toBe('navigate');
    expect(t.changeKind('drag', none, 'cell')).toBeUndefined();
    expect(changeKindOf('nothing')).toBeUndefined();
    expect(changeKindOf('view.fit', { 'view.fit': 'session' })).toBe('session');
    expect(Object.values(CommandChangeKinds).every((k) => ['patch', 'session', 'navigate'].includes(k))).toBe(true);
  });
});
