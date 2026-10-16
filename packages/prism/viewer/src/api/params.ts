/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

// The params the scene areas share: a `structure` and a `component` field name theirs by label or state ref (scope.ts
// resolves them with one rule), empty meaning the only one, or what the operation says. Each operation passes its own
// description where its meaning for "empty" differs.

/** A structure by label or state ref. Empty: the one in view, unless the operation's description says otherwise. */
export const StructureParam = (description = 'A structure by label or state ref; empty: the one in view') =>
  PD.Text('', { description });

/** A structure by label or state ref. Empty: every one (the operation acts on many). */
export const EveryStructureParam = () => StructureParam('A structure by label or state ref; empty: every one');

/** For an operation that finds a component among the structures' when none is named. */
export const InStructure = 'A structure by label or state ref; empty: every one (the component is found among them)';

/** A component by label or state ref. Empty: the only one, unless the operation's description says otherwise. */
export const ComponentParam = (description = 'A component by label or state ref; empty: the only one') =>
  PD.Text('', { description });

/**
 * A value an operation can leave as it is (operations.md §1): absent, it is kept; `checkValues` leaves it out (Mol*'s
 * `PD.Optional`). Its type stays the param's own, so read it in `run` as possibly undefined.
 */
export function OptionalParam<P extends PD.Any>(p: P): P {
  return { ...p, isOptional: true };
}

/** A switch's value (operations.md §1): off, or on with its settings, only those given set (the rest keep theirs). */
export type SwitchValue<S extends PD.Params = PD.Params> =
  | { readonly name: 'off'; readonly params?: Record<string, never> }
  | { readonly name: 'on'; readonly params?: Partial<PD.Values<S>> };

/** A switch's param, typed by its value (`SwitchValue`) for scripts; Mol*'s mapped on/off underneath. */
export type SwitchParamDef<S extends PD.Params> = Omit<PD.Mapped<any>, 'defaultValue'> & {
  readonly defaultValue: SwitchValue<S> | undefined;
};

/**
 * A switch with settings (operations.md §1), as Mol*'s own on/off props are: `{ name: 'off' }` or `{ name: 'on', params
 * }`, each setting optional (given: set, absent: keep); the whole switch absent keeps it as it is.
 */
export function SwitchParam<S extends PD.Params>(settings: S, info?: PD.Info): SwitchParamDef<S> {
  const optional: PD.Params = {};
  for (const [k, p] of Object.entries(settings)) optional[k] = PD.Optional(p) as PD.Any;
  return PD.Optional(
    PD.MappedStatic('on', { on: PD.Group(optional), off: PD.Group({}) }, info),
  ) as unknown as SwitchParamDef<S>;
}

/** The settings a switch turns on with: the current ones (or the defaults), the given ones over them. */
export function switchedOn<T extends object>(current: T, given: SwitchValue | undefined): T {
  return given?.name === 'on' && given.params ? { ...current, ...(given.params as Partial<T>) } : current;
}

/**
 * Params patched: a plain group merges one level (some of its values), anything else replaces; a value with a `kind`
 * (an iso value) replaces as a whole, so a relative one does not keep an absolute one's field.
 */
export function mergeParams(base: Record<string, any>, patch: Record<string, any>): Record<string, any> {
  const out: Record<string, any> = { ...base };
  for (const [k, v] of Object.entries(patch)) {
    out[k] = isPlain(v) && isPlain(out[k]) && !('kind' in v) ? { ...out[k], ...v } : v;
  }
  return out;
}

const isPlain = (x: unknown): x is Record<string, any> =>
  !!x && typeof x === 'object' && !Array.isArray(x) && Object.getPrototypeOf(x) === Object.prototype;
