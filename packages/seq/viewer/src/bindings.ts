/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The gesture bindings (_spec/README.md "Bindings"): `{ gesture, modifiers, target, action }`, one table per viewer.
// The action names a model method (continuous) or a command (discrete); the hit says the granularity. Pointer gestures
// only: keys go through `@molstar/ui`'s keymap. Built from a preset and overrides, checked against the actions that
// exist, conflicts reported when the table is built, never a race at run time.

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
  /** The modifiers held, exactly (none when absent); `'*'`: whatever is held, tried after the exact ones. */
  readonly modifiers?: readonly Modifier[] | '*';
  readonly target: BindingTarget;
  /** A method (`hover`, `pan`, `zoom`) or a command id; null in an override removes the binding. */
  readonly action: string | null;
}

/** The modifiers held, platform-aware. */
export interface ModifierState {
  shift?: boolean;
  mod?: boolean;
  alt?: boolean;
  space?: boolean;
}

/**
 * The modifiers of a pointer event: Cmd is `mod` on macOS, Ctrl elsewhere. With `out`, it is filled and returned
 * (one record reused for every move).
 */
export function modifiersOf(
  e: { shiftKey?: boolean; metaKey?: boolean; ctrlKey?: boolean; altKey?: boolean },
  mac: boolean,
  space = false,
  out: ModifierState = {},
): ModifierState {
  out.shift = !!e.shiftKey;
  out.mod = mac ? !!e.metaKey : !!e.ctrlKey;
  out.alt = !!e.altKey;
  out.space = space;
  return out;
}

const DataTargets: readonly HitTarget[] = ['cell', 'item'];
const ElementTargets: readonly HitTarget[] = ['cell', 'item', 'row', 'column', 'track'];

const each = (
  targets: readonly HitTarget[],
  gesture: PointerGesture,
  action: string,
  modifiers?: readonly Modifier[],
): GestureBinding[] => targets.map((target) => ({ gesture, target, action, ...(modifiers ? { modifiers } : {}) }));

/**
 * The selection's part of the grammar: added by the `select` behavior with the selection model and its commands; not
 * in the presets until they exist (a binding to nothing is reported).
 */
export const SelectionBindings: readonly GestureBinding[] = [
  ...each(ElementTargets, 'click', 'selection.set'),
  ...each(ElementTargets, 'click', 'selection.extend', ['shift']),
  ...each(ElementTargets, 'click', 'selection.toggle', ['mod']),
  { gesture: 'click', target: 'empty', action: 'selection.clear' },
  ...each(DataTargets, 'drag', 'drag-select'),
];

/**
 * The default grammar's hover and navigation. The hover holds with any modifier. Shift+wheel pans (sideways on a
 * vertical wheel); Alt is left unbound on purpose: the platform's and the host's (a context gesture, a copy-drag).
 */
export const DefaultBindings: readonly GestureBinding[] = [
  { gesture: 'hover', modifiers: '*', target: 'any', action: 'hover' },
  { gesture: 'drag', target: 'axis', action: 'pan' },
  { gesture: 'drag', target: 'empty', action: 'pan' },
  { gesture: 'drag', modifiers: ['space'], target: 'any', action: 'pan' },
  { gesture: 'wheel', target: 'any', action: 'pan' },
  { gesture: 'wheel', modifiers: ['shift'], target: 'any', action: 'pan' },
  { gesture: 'wheel', modifiers: ['mod'], target: 'any', action: 'zoom' },
  { gesture: 'pinch', target: 'any', action: 'zoom' },
  ...each(ElementTargets, 'double-click', 'view.focus'),
  { gesture: 'double-click', target: 'empty', action: 'view.fit' },
];

/** As genome browsers do: a drag on data pans too. */
export const GenomeBrowserBindings: readonly GestureBinding[] = [
  ...DefaultBindings,
  ...each(DataTargets, 'drag', 'pan'),
];

export const BindingPresets = {
  default: DefaultBindings,
  'genome-browser': GenomeBrowserBindings,
} as const satisfies Record<string, readonly GestureBinding[]>;

/** The actions the first slice's viewer has: the hover and navigation methods and the `navigate` commands. */
export const DefaultActions: readonly string[] = [
  'hover',
  'pan',
  'zoom',
  'view.fit',
  'view.focus',
  'view.zoom-in',
  'view.zoom-out',
  'view.pan-left',
  'view.pan-right',
];

export interface BindingTable {
  readonly bindings: readonly GestureBinding[];
  /** Conflicts, actions that do not exist, key gestures: found when the table is built. */
  readonly issues: readonly Issue[];
  /**
   * The action for a gesture with these modifiers on this kind of hit: the exact modifiers' binding (the target's,
   * else `any`'s), then the any-modifiers' one. Map lookups: no allocation.
   */
  resolve(gesture: PointerGesture, modifiers: ModifierState, target: HitTarget): string | undefined;
  /** The change kind of the action a gesture resolves to: undefined for a method or nothing. */
  changeKind(gesture: PointerGesture, modifiers: ModifierState, target: HitTarget): ChangeKind | undefined;
}

const ModifierOrder: readonly Modifier[] = ['shift', 'mod', 'alt', 'space'];
const AnyModifiers = -1;

function maskOf(mods: readonly Modifier[] | '*' | undefined): number {
  if (mods === '*') return AnyModifiers;
  let m = 0;
  for (const x of mods ?? []) m |= 1 << ModifierOrder.indexOf(x);
  return m;
}

const stateMask = (s: ModifierState) => (s.shift ? 1 : 0) | (s.mod ? 2 : 0) | (s.alt ? 4 : 0) | (s.space ? 8 : 0);

const keyOf = (gesture: string, mask: number, target: string) => `${gesture}|${mask}|${target}`;

export interface BindingTableOptions {
  /** The actions that exist: methods and command ids (`DefaultActions`). A binding to another is reported and left out. */
  readonly actions?: readonly string[];
  /** Change kinds of commands the kind table lacks (an extension's). */
  readonly kinds?: Readonly<Record<string, ChangeKind>>;
}

/**
 * A table from a preset and overrides (an application's, then a person's), each layer replacing the bindings of the
 * same gesture, modifiers and target before it (`action: null` removes one). Two bindings of one key within a layer
 * conflict; an action that does not exist is reported and left out; key gestures are reported and left to the keymap.
 */
export function bindingTable(
  preset: readonly GestureBinding[] = DefaultBindings,
  overrides: readonly (readonly GestureBinding[])[] = [],
  options: BindingTableOptions = {},
): BindingTable {
  const actions = new Set(options.actions ?? DefaultActions);
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
      if (!actions.has(b.action)) {
        issues.push(
          issueOf(ViewerIssues, 'seq.viewer.binding.action', `No action ${b.action} for ${describe(b)}`, {
            hint: ViewerIssues.codes['seq.viewer.binding.action'].hint,
          }),
        );
        table.delete(key);
        continue;
      }
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
    const byMask = index.get(gesture);
    if (!byMask) return undefined;
    const exact = byMask.get(stateMask(modifiers));
    const hit = exact?.get(target) ?? exact?.get('any');
    if (hit !== undefined) return hit;
    const any = byMask.get(AnyModifiers);
    return any?.get(target) ?? any?.get('any');
  };
  return {
    bindings: [...table.values()],
    issues,
    resolve,
    changeKind: (gesture, modifiers, target) => {
      const action = resolve(gesture, modifiers, target);
      return action === undefined || isMethodAction(action) ? undefined : changeKindOf(action, options.kinds);
    },
  };
}

function describe(b: GestureBinding): string {
  const mods = b.modifiers === '*' ? '(any)+' : b.modifiers?.length ? `${b.modifiers.join('+')}+` : '';
  return `${mods}${b.gesture} on ${b.target}`;
}
