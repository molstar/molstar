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

/** A component by label or state ref. Empty: the only one, unless the operation's description says otherwise. */
export const ComponentParam = (description = 'A component by label or state ref; empty: the only one') =>
  PD.Text('', { description });

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
