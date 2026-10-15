/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The behaviors' bindings for a help panel (seqspec:viewer#behaviors): every binding param of a built-in behavior is
// described, and `ctx.behaviors.bindings()` lists them from the current params as Mol*'s `Binding`, which Mol*'s
// `Binding.format` renders.

import { Binding as MolBinding } from '@molstar/core/util/binding';
import { describe, expect, it } from 'vitest';
import { case3 } from '../../core/_test/format/cases.js';
import { SeqViewerBehaviors, type SeqBehaviorProvider } from '../src/index.js';
import { driven } from './input-fixtures.js';

/** The params' paths whose value is a list of triggers (a pointer's `button` or modifiers, a key's `key` or `code`). */
function bindingPaths(params: object, prefix = ''): string[] {
  const out: string[] = [];
  for (const [k, v] of Object.entries(params)) {
    const path = prefix ? `${prefix}.${k}` : k;
    const trigger = (x: unknown) =>
      x !== null &&
      typeof x === 'object' &&
      Object.keys(x).every((f) => /^(button|key|code|shift|ctrl|alt|meta)$/.test(f));
    if (Array.isArray(v) && v.length > 0 && v.every(trigger)) out.push(path);
    else if (v !== null && typeof v === 'object' && !Array.isArray(v)) out.push(...bindingPaths(v, path));
  }
  return out;
}

describe('bindings', () => {
  it('every binding param of a built-in behavior is described', () => {
    const missing: string[] = [];
    for (const b of Object.values(SeqViewerBehaviors) as SeqBehaviorProvider[])
      for (const path of bindingPaths(b.params)) if (!b.bindings?.[path]) missing.push(`${b.name}.${path}`);
    expect(missing).toEqual([]);
  });

  it("lists them from the current params as Mol*'s Binding, which Mol* formats", async () => {
    const t = await driven(case3());
    t.ctx.behaviors.get('select')?.setParams({ toggle: [{ button: 'primary', meta: true }] });
    const all = t.ctx.behaviors.bindings();
    const toggle = all.find((b) => b.behavior === 'select' && b.param === 'toggle');
    expect(toggle?.binding).toEqual({
      triggers: [{ buttons: 1, modifiers: { shift: false, control: false, alt: false, meta: true } }],
      action: 'toggle',
      description: 'Add or remove the element with ${triggers}',
    });
    const menu = all.find((b) => b.behavior === 'keyboard' && b.param === 'bindings.menu');
    expect(menu?.binding.triggers).toEqual([
      { code: 'F10', modifiers: { shift: true, control: false, alt: false, meta: false } },
      { code: 'ContextMenu', modifiers: { shift: false, control: false, alt: false, meta: false } },
    ]);
    expect(MolBinding.format(menu?.binding as MolBinding)).toBe(
      "Ask for the track's menu with <i>f10 + shift key or context menu</i>",
    );
    const zoom = all.find((b) => b.behavior === 'wheel' && b.param === 'zoom');
    expect(zoom?.binding.triggers[0]).toEqual({
      buttons: undefined,
      modifiers: { shift: false, control: true, alt: false, meta: false },
    });
    // An unset trigger (`dragPan`'s `dominantAxis`) is not listed.
    expect(all.some((b) => b.param === 'dominantAxis')).toBe(false);
    t.ctx.behaviors.get('keyboard')?.remove();
    expect(t.ctx.behaviors.bindings().some((b) => b.behavior === 'keyboard')).toBe(false);
    t.ctx.dispose();
  });
});
