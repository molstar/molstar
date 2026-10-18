/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { Binding } from '@molstar/core/util/binding';
import { DefaultTrackballBindings } from '@molstar/graphics/canvas3d/controls/trackball';
import { CameraControls } from '@molstar/plugin/behavior/dynamic/camera';
import { SnapshotControls } from '@molstar/plugin/behavior/dynamic/state';
import { describe, expect, it } from 'vitest';
import { CanvasBindings } from '../bindings.js';

// Help's keys are written out (`ui/bindings.ts`): every one must be a key the plugin binds by default, and every key
// the plugin binds (but the gamepad's) must be listed.

const Arrows: Record<string, string> = { ArrowUp: '↑', ArrowDown: '↓', ArrowLeft: '←', ArrowRight: '→' };

/** A key trigger as Help writes it: the modifiers, then the key ('⇧', '↑'). */
function keysOf(t: Binding.Trigger): string[] | undefined {
  const code = t.code;
  if (!code || code.startsWith('Gamepad')) return undefined;
  const m = t.modifiers;
  const mods = [m?.shift ? '⇧' : '', m?.control ? '⌃' : '', m?.alt ? '⌥' : '', m?.meta ? '⌘' : ''].filter(Boolean);
  const key = Arrows[code] ?? (code === 'ShiftLeft' ? '⇧' : code.replace(/^Key/, ''));
  return [...mods, key];
}

function defaultsOf(behavior: { definition: { params?: (a: any, ctx: any) => PD.Params } }) {
  return PD.getDefaultValues(behavior.definition.params!(undefined, undefined)).bindings as Record<string, Binding>;
}

describe('the Help bindings', () => {
  const plugin = [DefaultTrackballBindings, defaultsOf(CameraControls), defaultsOf(SnapshotControls)].flatMap((b) =>
    Object.values(b as Record<string, Binding>).flatMap((binding) =>
      binding.triggers.flatMap((t) => {
        const keys = keysOf(t);
        return keys ? [keys.join(' ')] : [];
      }),
    ),
  );
  const listed = CanvasBindings.flatMap((g) => g.bindings)
    .filter((b) => b.input.every((k) => k.length <= 5 && !/drag|click|scroll|hover/i.test(k)))
    .map((b) => b.input.join(' '));

  it('lists every key the plugin binds', () => {
    expect(plugin.length).toBeGreaterThan(15);
    expect(plugin.filter((k) => !listed.includes(k))).toEqual([]);
  });

  it('lists no key the plugin does not bind', () => {
    expect(listed.filter((k) => !plugin.includes(k))).toEqual([]);
  });
});
