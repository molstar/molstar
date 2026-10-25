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
  lookCellsUnder,
  cellTakes,
  cellLookValue,
  structureReprsOf,
} from '../api/look.js';
import { componentName, isFailure, labelOf } from '../api/scope.js';
import type { LookOverrideSummary, NodeLookModel, SceneLookMark, SceneLookModel } from '../model/look.js';
import type { SourceScope } from './state.js';

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
 * sent, so a value another surface set meanwhile is not written back.
 */
function batched(scope: SourceScope | undefined, run: (changed: Record<string, unknown>) => void, ms = 150) {
  let timer: ReturnType<typeof setTimeout> | undefined;
  let pending: Record<string, unknown> = {};
  scope?.add(() => clearTimeout(timer));
  return (values: Record<string, unknown>, name: string) => {
    pending[name] = values[name];
    clearTimeout(timer);
    timer = setTimeout(() => {
      const changed = pending;
      pending = {};
      run(changed);
    }, ms);
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
      list.push(reprName(plugin, r));
    }
  }
  return LookKeyNames.filter((k) => names.has(k)).map((k) => ({ key: k, label: labelOfKey(k), names: names.get(k)! }));
}

function reprName(plugin: PluginContext, r: StructureRepresentationRef): string {
  const type = (r.cell.transform.params as { type: { name: string } }).type.name;
  const name = plugin.representation.structure.registry.get(type)?.label ?? type;
  return `${name} of ${componentName(r.component)} (${labelOf({ plugin }, r.component.structure)})`;
}

const sameSummary = (a: readonly LookOverrideSummary[], b: readonly LookOverrideSummary[]) =>
  a.length === b.length &&
  a.every(
    (x, i) =>
      x.key === b[i].key && x.names.length === b[i].names.length && x.names.every((n, j) => n === b[i].names[j]),
  );

/** The scene defaults (the Structure panel's sliders): `look.set` of the scene, the overrides kept, Apply to all. */
export function createSceneLookModel(scope: SourceScope, plugin: PluginContext, api: ViewerApi): SceneLookModel {
  const params = sceneFormParams();
  const form = new ParamsFormModel(params, sceneValues(plugin) as never, {
    onChange: batched(scope, (changed) => void api.run('look.set', changed)) as never,
  });
  const overrides = observableState<readonly LookOverrideSummary[]>([], { compare: sameSummary });
  // Read once per burst of changes (a drag's updates, a preset's many cells), after it.
  let queued = false;
  const read = () => {
    queued = false;
    overrides.set(summarize(plugin, everyRepresentation({ plugin })));
    const values = sceneValues(plugin);
    if (!PD.areEqual(params, form.values.value, values)) form.load(values as never);
  };
  const queue = () => {
    if (queued) return;
    queued = true;
    queueMicrotask(read);
  };
  scope.add(plugin.managers.structure.component.events.optionsUpdated.subscribe(queue));
  scope.add(plugin.state.data.events.changed.subscribe(queue));
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
  const cells = () => {
    const r = lookCellsUnder({ plugin }, ref);
    return isFailure(r) ? [] : r;
  };
  const found = cells();
  if (!found.length) return undefined;
  const all = lookFormParams() as Record<string, PD.Any>;
  const params: Record<string, PD.Any> = {};
  for (const k of LookKeyNames) if (all[k] && found.some((c) => cellTakes(c, k))) params[k] = all[k];
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
  const form = new ParamsFormModel(params, values() as never, {
    onChange: batched(undefined, (changed) => void api.run('look.set', { node: ref, ...changed })) as never,
  });
  const reprs = () => structureReprsOf({ plugin }, cells());
  const own = observableState<readonly LookOverrideSummary[]>(summarize(plugin, reprs()), { compare: sameSummary });
  return {
    form,
    own: own.asReadonly,
    reset: (key) => void api.run('look.reset', { node: ref, keys: key ? [key] : [] }),
    refresh: () => {
      own.set(summarize(plugin, reprs()));
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
