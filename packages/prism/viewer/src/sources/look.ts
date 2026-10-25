/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { StructureComponentManager } from '@molstar/plugin/state/manager/structure/component';
import type { StructureRepresentationRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { observableState } from '@molstar/prism-ui/core';
import { ParamsFormModel } from '@molstar/prism-ui/forms/params-form';
import type { ViewerApi } from '../api/area.js';
import {
  everyRepresentation,
  type LookKey,
  LookKeyNames,
  LookLabels,
  lookDefaults,
  ownLookValues,
  reprLabel,
  lookCellsUnder,
  cellTakes,
  cellLookValue,
  structureReprsOf,
} from '../api/look.js';
import { isFailure } from '../api/scope.js';
import type { LookOverrideSummary, NodeLookModel, SceneLookMark, SceneLookModel } from '../model/look.js';
import { onceATask, sameItems, type SourceScope } from './state.js';
import { merge } from 'rxjs';
import type { StateObjectCell } from '@molstar/core/state';

// The look's view models (api.md §2, "Look"): the scene defaults as the Structure panel's sliders show them, a node's
// own values as its inspector's sliders show them, and the tree rows' marks. Each change is a `look.*` run.

const P = StructureComponentManager.OptionsParams;

/** A look value's name for people: "Ignore light". */
const labelOfKey = (k: string) => {
  const text = LookLabels[k as LookKey] ?? k;
  return text[0].toUpperCase() + text.slice(1);
};

/** The look's params by `look.set`'s names, labeled as the panel shows them. */
function lookFormParams() {
  return {
    hydrogens: { ...P.hydrogens, label: 'Hydrogens' },
    quality: { ...P.visualQuality, label: 'Visual Quality' },
    ignoreLight: { ...P.ignoreLight, label: 'Ignore Light' },
    material: { ...P.materialStyle, label: 'Material Style' },
    solidSurface: { ...P.solidSurface, label: 'Solid Surface' },
    interior: P.interior,
  };
}

/**
 * The scene defaults' form: the look but the solid surface (the Clipping section's Settings), and the interactions'
 * settings. The clip objects are the Clipping section's cuts, the animation the Animation card's Wiggle.
 */
function sceneFormParams() {
  const { solidSurface: _, ...look } = lookFormParams();
  return { ...look, interactions: P.interactions };
}

/** The component manager's options by the form's names. */
function sceneValues(plugin: PluginContext) {
  const { hydrogens, quality, ignoreLight, material, interior } = lookDefaults({ plugin });
  const interactions = plugin.managers.structure.component.state.options.interactions;
  return { hydrogens, quality, ignoreLight, material, interior, interactions };
}

/**
 * Collects changed fields while a drag goes on, then runs once (150 ms after the last change): only what changed is
 * sent, so a value another surface set meanwhile is not written back. `busy` while a change waits or runs: the form is
 * not loaded from the scene then (it would snap back to the value before the drag). The timer goes with `scope`.
 */
function batched(
  scope: SourceScope | undefined,
  run: (changed: Record<string, unknown>) => Promise<unknown>,
  ms = 150,
) {
  let timer: ReturnType<typeof setTimeout> | undefined;
  let pending: Record<string, unknown> | undefined;
  let running = 0;
  scope?.add(() => clearTimeout(timer));
  return {
    onChange: (values: Record<string, unknown>, name: string) => {
      (pending ??= {})[name] = values[name];
      clearTimeout(timer);
      timer = setTimeout(() => {
        const changed = pending!;
        pending = undefined;
        running++;
        void run(changed).finally(() => running--);
      }, ms);
    },
    busy: () => pending !== undefined || running > 0,
    cancel: () => {
      clearTimeout(timer);
      pending = undefined;
    },
  };
}

/**
 * A representation's own values (their keys), cached by its transform and the options: a state change that leaves both
 * (a trajectory's frame, another cell's update) reads nothing again.
 */
const ownCache = new WeakMap<object, { options: object; keys: readonly LookKey[] }>();
function ownKeys(plugin: PluginContext, r: StructureRepresentationRef): readonly LookKey[] {
  const options = plugin.managers.structure.component.state.options;
  const hit = ownCache.get(r.cell.transform);
  if (hit && hit.options === options) return hit.keys;
  const keys = Object.keys(ownLookValues({ plugin }, r)) as LookKey[];
  ownCache.set(r.cell.transform, { options, keys });
  return keys;
}

/** The overrides among these representations, by value. */
function summarize(plugin: PluginContext, reprs: readonly StructureRepresentationRef[]): LookOverrideSummary[] {
  const names = new Map<LookKey, string[]>();
  for (const r of reprs) {
    for (const k of ownKeys(plugin, r)) {
      let list = names.get(k);
      if (!list) names.set(k, (list = []));
      list.push(reprLabel({ plugin }, r));
    }
  }
  return LookKeyNames.filter((k) => names.has(k)).map((k) => ({ key: k, label: labelOfKey(k), names: names.get(k)! }));
}

const sameSummary = (a: readonly LookOverrideSummary[], b: readonly LookOverrideSummary[]) =>
  sameItems(a, b, (x, y) => x.key === y.key && sameItems(x.names, y.names));

/**
 * The scene defaults (the Structure panel's sliders): `look.set` of the scene, the overrides kept, Apply to all. Read
 * again once a burst of what can change them: the options, the structures' hierarchy (a representation's params, one
 * added or removed); not every state change (a trajectory's frames). The form follows unless a change of its own waits.
 */
export function createSceneLookModel(scope: SourceScope, plugin: PluginContext, api: ViewerApi): SceneLookModel {
  const params = sceneFormParams();
  const changes = batched(scope, (changed) => api.run('look.set', changed));
  const form = new ParamsFormModel(params, sceneValues(plugin) as never, { onChange: changes.onChange as never });
  const overrides = observableState<readonly LookOverrideSummary[]>([], { compare: sameSummary });
  const read = () => {
    overrides.set(summarize(plugin, everyRepresentation({ plugin })));
    if (changes.busy()) return;
    const values = sceneValues(plugin);
    if (!PD.areEqual(params, form.values.value, values)) form.load(values as never);
  };
  scope.add(
    onceATask(
      merge(
        plugin.managers.structure.component.events.optionsUpdated,
        plugin.managers.structure.hierarchy.behaviors.selection,
      ),
    ).subscribe(read),
  );
  read();
  return {
    form,
    overrides: overrides.asReadonly,
    applyToAll: (key) => void api.run('look.reset', { keys: [key] }),
  };
}

/**
 * A node's look (its inspector's sliders): the values of the representations under it (structures', volumes',
 * particles', shapes'), a change `look.set` of the node, Use default. Only the values some of them take are in the form
 * (hydrogens: a structure's); only structure representations follow the scene defaults, so only they have values of
 * their own.
 */
export function createNodeLookModel(plugin: PluginContext, api: ViewerApi, ref: string): NodeLookModel | undefined {
  // The cells under the node, found again only when the tree's shape changed (not at a param update: a frame)
  let found: { children: unknown; cells: StateObjectCell[] } | undefined;
  const cells = () => {
    const children = plugin.state.data.tree.children;
    if (found?.children !== children) {
      const r = lookCellsUnder({ plugin }, ref);
      found = { children, cells: isFailure(r) ? [] : r };
    }
    return found.cells;
  };
  if (!cells().length) return undefined;
  const all = lookFormParams() as Record<string, PD.Any>;
  const params: Record<string, PD.Any> = {};
  for (const k of LookKeyNames) if (all[k] && cells().some((c) => cellTakes(c, k))) params[k] = all[k];
  // each value as the first representation that takes it has it
  const values = () => {
    const now = cells();
    const out = lookDefaults({ plugin }) as Record<string, unknown>;
    for (const k of Object.keys(params) as LookKey[]) {
      const cell = now.find((c) => cellTakes(c, k));
      const v = cell ? cellLookValue(cell, k) : undefined;
      if (v !== undefined) out[k] = v;
    }
    return out;
  };
  // a change waiting when the node is gone is dropped (the inspector closed with it)
  const changes = batched(undefined, (changed) =>
    plugin.state.data.cells.has(ref) ? api.run('look.set', { node: ref, ...changed }) : Promise.resolve(),
  );
  const form = new ParamsFormModel(params, values() as never, { onChange: changes.onChange as never });
  const reprs = () => structureReprsOf({ plugin }, cells());
  const own = observableState<readonly LookOverrideSummary[]>(summarize(plugin, reprs()), { compare: sameSummary });
  // What the last refresh read: the options and the cells' transforms (the same: nothing to read again)
  let seen: readonly unknown[] = [];
  return {
    form,
    own: own.asReadonly,
    reset: (key) => void api.run('look.reset', { node: ref, keys: key ? [key] : [] }),
    refresh: () => {
      const now = [plugin.managers.structure.component.state.options, ...cells().map((c) => c.transform)];
      if (sameItems(now, seen)) return;
      seen = now;
      own.set(summarize(plugin, reprs()));
      if (changes.busy()) return;
      const next = values();
      if (!PD.areEqual(params, form.values.value, next)) form.load(next as never);
    },
  };
}

/** A representation row's mark: the values it keeps as its own. */
export function representationLookMark(
  plugin: PluginContext,
  r: StructureRepresentationRef,
): SceneLookMark | undefined {
  const keys = ownKeys(plugin, r);
  if (!keys.length) return undefined;
  const list = keys.map((k) => LookLabels[k]);
  const text = list.length === 1 ? list[0] : `${list.slice(0, -1).join(', ')} and ${list[list.length - 1]}`;
  return { kind: 'own', text: `Its own ${text}, not the scene default` };
}

/** A structure's or component's row mark: how many representations under it have looks of their own. */
export function groupLookMark(
  plugin: PluginContext,
  reprs: readonly StructureRepresentationRef[],
): SceneLookMark | undefined {
  let n = 0;
  for (const r of reprs) if (ownKeys(plugin, r).length) n++;
  if (!n) return undefined;
  return {
    kind: 'inside',
    text: n === 1 ? 'A representation here has its own look' : `${n} representations here have their own look`,
  };
}

/** Whether two marks are the same (a row writes only a change). */
export const sameLookMark = (a: SceneLookMark | undefined, b: SceneLookMark | undefined) =>
  a === b || (!!a && !!b && a.kind === b.kind && a.text === b.text);
