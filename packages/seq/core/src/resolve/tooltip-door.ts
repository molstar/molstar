/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type * as Tables from './tooltip-tables.js';

// The tooltips' door (_spec/README.md "Resolve", Tooltip tables, Loaded when used): a tooltip's resolve past its
// declaration (its tables as data, here), and the parts the other layers register (`withTooltips`: the view's pass,
// render's lookup and wants, a viewer's card provider) load the first time a spec holds a `tooltip`, a model waiting
// for them as for the genome chunk (`StaticModel.loading`), so a spec's first hover finds them landed. Why: only some
// specs carry tooltips; a reader of one without pays nothing for them.

/** Every registered part landed (a part registered later sets it back: the next model with a tooltip loads it). */
let landed = false;
let loading: Promise<void> | undefined;
const Parts: (() => Promise<unknown>)[] = [];

/** The tables' module when it has landed (`tables`), else none: resolve reads no tooltip table until then. */
export const TooltipTables: { tables?: typeof Tables } = {};

/** A part that loads with the tooltips (a layer's hooks set by its `load`); a later one loads with the next model. */
export function withTooltips(load: () => Promise<unknown>): void {
  Parts.push(load);
  landed = false;
}

/**
 * Loads the tooltips' parts, once while in flight. Never rejects: a failed import leaves them unloaded, and the next
 * model holding a tooltip tries again (its tooltips show nothing meanwhile, never wrong).
 */
export function loadTooltips(): Promise<void> {
  if (landed) return Promise.resolve();
  const n = Parts.length;
  loading ??= Promise.all([import('./tooltip-tables.js'), ...Parts.map((load) => load())]).then(
    ([m]) => {
      TooltipTables.tables = m;
      landed = Parts.length === n;
      loading = undefined;
    },
    () => {
      loading = undefined;
    },
  );
  return loading;
}

/** Whether a spec holds a tooltip (under `views`, or a view). */
export function hasTooltip(spec: Node): boolean {
  for (const c of spec.children ?? []) {
    if (c.kind !== 'views') continue;
    for (const n of c.children ?? []) {
      if (n.kind === 'tooltip') return true;
      if (n.kind === 'view') for (const t of n.children ?? []) if (t.kind === 'tooltip') return true;
    }
  }
  return false;
}

/** The load a model holding a tooltip waits for, while the parts have not landed; undefined when nothing is to load. */
export function tooltipsLoading(spec: Node): Promise<void> | undefined {
  return landed || !hasTooltip(spec) ? undefined : loadTooltips();
}
