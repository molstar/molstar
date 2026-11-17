/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The gesture bindings (_spec/README.md "Bindings"): `{ gesture, modifiers, target, action }`, one
// table per viewer. The action names a model method (continuous) or a command (discrete); the hit says the
// granularity. Pointer gestures only: keys go through `@molstar/ui`'s keymap. Built from a preset and overrides,
// conflicts reported when the table is built, never a race at run time.

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { type ChangeKind, changeKindOf, isMethodAction } from './commands.js';
import { ViewerIssues } from './issues.js';
import type { HitTarget } from './models/hover.js';

export type PointerGesture = 'hover' | 'click' | 'double-click' | 'drag' | 'wheel' | 'pinch';
export type Gesture = PointerGesture | 'key';

/** `mod` is Cmd on macOS, Ctrl elsewhere (`modifiersOf`); `space` is the key held during a drag. */
export type Modifier = 'shift' | 'mod' | 'alt' | 'space';

/** A hit's kind, or any: a binding for a kind beats one for any. */
export type BindingTarget = HitTarget | 'any';

export interface GestureBinding {
  readonly gesture: Gesture;
  readonly modifiers?: readonly Modifier[];
  readonly target: BindingTarget;
  /** A method (`hover`, `pan`, `zoom`, `drag-select`) or a command id; null in an override removes the binding. */
  readonly action: string | null;
}

/** The modifiers held, platform-aware. */
export interface ModifierState {
  readonly shift?: boolean;
  readonly mod?: boolean;
  readonly alt?: boolean;
  readonly space?: boolean;
}

/** The modifiers of a pointer event: Cmd is `mod` on macOS, Ctrl elsewhere. */
export function modifiersOf(
  e: { shiftKey?: boolean; metaKey?: boolean; ctrlKey?: boolean; altKey?: boolean },
  mac: boolean,
  space = false,
): ModifierState {
  return { shift: !!e.shiftKey, mod: mac ? !!e.metaKey : !!e.ctrlKey, alt: !!e.altKey, space };
}

const DataTargets: readonly HitTarget[] = ['cell', 'item'];
const ElementTargets: readonly HitTarget[] = ['cell', 'item', 'row', 'column', 'track'];

const each = (
  targets: readonly HitTarget[],
  gesture: PointerGesture,
  action: string,
  modifiers?: readonly Modifier[],
): GestureBinding[] => targets.map((target) => ({ gesture, target, action, ...(modifiers ? { modifiers } : {}) }));

/** The default grammar. */
export const DefaultBindings: readonly GestureBinding[] = [
  { gesture: 'hover', target: 'any', action: 'hover' },
  ...each(ElementTargets, 'click', 'selection.set'),
  ...each(ElementTargets, 'click', 'selection.extend', ['shift']),
  ...each(ElementTargets, 'click', 'selection.toggle', ['mod']),
  { gesture: 'click', target: 'empty', action: 'selection.clear' },
  ...each(DataTargets, 'drag', 'drag-select'),
  { gesture: 'drag', target: 'axis', action: 'pan' },
  { gesture: 'drag', target: 'empty', action: 'pan' },
  { gesture: 'drag', modifiers: ['space'], target: 'any', action: 'pan' },
  { gesture: 'wheel', target: 'any', action: 'pan' },
  { gesture: 'wheel', modifiers: ['mod'], target: 'any', action: 'zoom' },
  { gesture: 'pinch', target: 'any', action: 'zoom' },
  ...each(ElementTargets, 'double-click', 'view.focus'),
  { gesture: 'double-click', target: 'empty', action: 'view.fit' },
];

/** As genome browsers do: a drag on data pans; Shift+drag selects. */
export const GenomeBrowserBindings: readonly GestureBinding[] = [
  ...DefaultBindings.filter((b) => !(b.gesture === 'drag' && b.action === 'drag-select')),
  ...each(DataTargets, 'drag', 'pan'),
  ...each(DataTargets, 'drag', 'drag-select', ['shift']),
];

export const BindingPresets = {
  default: DefaultBindings,
  'genome-browser': GenomeBrowserBindings,
} as const satisfies Record<string, readonly GestureBinding[]>;

export interface BindingTable {
  readonly bindings: readonly GestureBinding[];
  /** Conflicts, unknown actions, key gestures: found when the table is built. */
  readonly issues: readonly Issue[];
  /** The action for a gesture with these modifiers on this kind of hit: the target's binding, else `any`'s. */
  resolve(gesture: PointerGesture, modifiers: ModifierState, target: HitTarget): string | undefined;
  /** The change kind of the action a gesture resolves to: undefined for a method or nothing. */
  changeKind(gesture: PointerGesture, modifiers: ModifierState, target: HitTarget): ChangeKind | undefined;
}

const ModifierOrder: readonly Modifier[] = ['shift', 'mod', 'alt', 'space'];

function maskOf(mods: readonly Modifier[] | ModifierState | undefined): number {
  if (!mods) return 0;
  let m = 0;
  if (Array.isArray(mods)) {
    for (const x of mods as readonly Modifier[]) m |= 1 << ModifierOrder.indexOf(x);
    return m;
  }
  const s = mods as ModifierState;
  return (s.shift ? 1 : 0) | (s.mod ? 2 : 0) | (s.alt ? 4 : 0) | (s.space ? 8 : 0);
}

const keyOf = (gesture: string, mask: number, target: string) => `${gesture}|${mask}|${target}`;

/**
 * A table from a preset and overrides (an application's, then a person's), each layer replacing the bindings of the
 * same gesture, modifiers and target before it (`action: null` removes one). Two bindings of one key within a layer
 * conflict; a name that is no method and no command known to `commands` (the change-kind table and `extra`) is
 * reported; key gestures are reported and left to the keymap.
 */
export function bindingTable(
  preset: readonly GestureBinding[] = DefaultBindings,
  overrides: readonly (readonly GestureBinding[])[] = [],
  options: { commands?: Readonly<Record<string, ChangeKind>> } = {},
): BindingTable {
  const issues: Issue[] = [];
  const table = new Map<string, GestureBinding>();
  for (const layer of [preset, ...overrides]) {
    const seen = new Set<string>();
    for (const b of layer) {
      if (b.gesture === 'key') {
        issues.push(issueOf(ViewerIssues, 'seq.viewer.binding.key', `A key binding for ${b.action}`));
        continue;
      }
      const key = keyOf(b.gesture, maskOf(b.modifiers), b.target);
      if (seen.has(key))
        issues.push(issueOf(ViewerIssues, 'seq.viewer.binding.conflict', `Two bindings for ${describe(b)}`));
      seen.add(key);
      if (b.action === null) {
        table.delete(key);
        continue;
      }
      if (!isMethodAction(b.action) && changeKindOf(b.action, options.commands) === undefined)
        issues.push(
          issueOf(ViewerIssues, 'seq.viewer.binding.action', `Unknown action ${b.action} for ${describe(b)}`),
        );
      table.set(key, b);
    }
  }
  // a lookup per gesture, mask and target, so resolving allocates nothing
  const index = new Map<string, Map<number, Map<string, string>>>();
  for (const b of table.values()) {
    let byMask = index.get(b.gesture);
    if (!byMask) index.set(b.gesture, (byMask = new Map()));
    const mask = maskOf(b.modifiers);
    let byTarget = byMask.get(mask);
    if (!byTarget) byMask.set(mask, (byTarget = new Map()));
    byTarget.set(b.target, b.action!);
  }
  const resolve = (gesture: PointerGesture, modifiers: ModifierState, target: HitTarget) => {
    const byTarget = index.get(gesture)?.get(maskOf(modifiers));
    return byTarget?.get(target) ?? byTarget?.get('any');
  };
  return {
    bindings: [...table.values()],
    issues,
    resolve,
    changeKind: (gesture, modifiers, target) => {
      const action = resolve(gesture, modifiers, target);
      return action === undefined || isMethodAction(action) ? undefined : changeKindOf(action, options.commands);
    },
  };
}

function describe(b: GestureBinding): string {
  const mods = b.modifiers?.length ? `${b.modifiers.join('+')}+` : '';
  return `${mods}${b.gesture} on ${b.target}`;
}
