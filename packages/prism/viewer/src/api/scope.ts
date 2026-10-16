/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransform } from '@molstar/core/state';
import type { StructureComponentParams } from '@molstar/plugin/state/helpers/structure-component';
import type { StructureComponentRef, StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { entryOf, structureNames, workingStructure } from '../sources/naming.js';
import { type ApiContext, type ApiFailure, ApiResult } from './area.js';

// Which things an operation acts on (targets.md §1, "the structure"): the one rule that turns a name into one of a
// list, for structures, components, volumes and particle lists alike, and the structures an operation acts on by it.
// A name is a state ref (exact), else a name for people (case-insensitive); a name that matches several is `ambiguous`,
// with the candidates (name and state ref) in the hint; no name is the only one there is, else `ambiguous`; none is
// `not-found`, listing what exists. Shared by the scene areas.

/** One failure type for the areas: `ApiResult.isFailure`. */
export const isFailure = ApiResult.isFailure;

// --- the rule ---------------------------------------------------------------------------------------------------

/** How things of one kind are told apart by name. */
export interface NameRule<T> {
  /** What they are, for messages: "structure", "component", "particle list". */
  readonly what: string;
  /** The identifiers that match exactly: the state ref, and a key where there is one. */
  ids(x: T): readonly string[];
  /** The names for people, most specific first: the first position with a match decides (case-insensitive). */
  names(x: T): readonly (string | undefined)[];
  /** A candidate in a hint: its name and state ref. */
  describe(x: T): string;
  /** When there are none at all. */
  readonly none: { readonly message: string; readonly hint?: string };
  /** The hint's lead-in to what exists on a name that is not found. */
  readonly listed?: string;
  /** Said after "Which X? n" when no name is given and there are several. */
  readonly several?: string;
  /** Said after "Give its state ref" when a name matches several. */
  readonly more?: string;
}

/** One of `list` by `name` (see above); `undefined`/empty: the only one, else `ambiguous`. */
export function pickNamed<T>(list: readonly T[], name: string | undefined, rule: NameRule<T>): T | ApiFailure {
  if (!list.length) return ApiResult.fail('not-found', rule.none.message, rule.none.hint);
  const all = () => list.map(rule.describe).join(', ');
  if (!name) {
    if (list.length === 1) return list[0];
    return ApiResult.fail(
      'ambiguous',
      `Which ${rule.what}? ${list.length} ${rule.several ?? 'are loaded'}.`,
      `Name one: ${all()}.`,
    );
  }
  let found = list.filter((x) => rule.ids(x).includes(name));
  if (!found.length) {
    const lower = name.trim().toLowerCase();
    const names = list.map(rule.names);
    const tiers = names.reduce((n, x) => Math.max(n, x.length), 0);
    for (let i = 0; i < tiers && !found.length; i++) {
      found = list.filter((_, j) => names[j][i]?.toLowerCase() === lower);
    }
  }
  if (found.length === 1) return found[0];
  if (found.length > 1) {
    return ApiResult.fail(
      'ambiguous',
      `Which '${name}'? ${found.length} ${rule.what}s have that name.`,
      `Give its state ref${rule.more ?? ''}: ${found.map(rule.describe).join(', ')}.`,
    );
  }
  return ApiResult.fail('not-found', `There is no ${rule.what} '${name}'.`, `${rule.listed ?? 'Loaded'}: ${all()}.`);
}

// --- structures -------------------------------------------------------------------------------------------------

export function loadedStructures(ctx: Pick<ApiContext, 'plugin'>): readonly StructureRef[] {
  return ctx.plugin.managers.structure.hierarchy.current.structures;
}

/** A structure as every read names it (api.md §1 "Words"): its state ref and its name for people. */
export interface StructureId {
  readonly ref: string;
  readonly label: string;
}

export function structureId(ctx: Pick<ApiContext, 'plugin'>, s: StructureRef): StructureId {
  return { ref: s.cell.transform.ref, label: labelOf(ctx, s) };
}

/** The structure's name for people (`sources/naming.ts`): its entry, told apart when two share one. */
export function labelOf(ctx: Pick<ApiContext, 'plugin'>, s: StructureRef): string {
  return structureNames(loadedStructures(ctx)).get(s.cell.transform.ref) ?? s.cell.transform.ref;
}

function structureRule(list: readonly StructureRef[]): NameRule<StructureRef> {
  const names = structureNames(list);
  const labelOfRef = (s: StructureRef) => names.get(s.cell.transform.ref) ?? s.cell.transform.ref;
  return {
    what: 'structure',
    ids: (s) => [s.cell.transform.ref],
    // The name told apart, the entry, the plugin's label ("Assembly 1").
    names: (s) => [names.get(s.cell.transform.ref), entryOf(s), s.cell.obj?.label],
    describe: (s) => `${labelOfRef(s)} (${s.cell.transform.ref})`,
    none: { message: 'Nothing is loaded.', hint: 'Load a structure first.' },
  };
}

/**
 * One structure (targets.md §1, "the structure"): the named one, else the one being worked with (focused, else the one
 * holding the selection), else the only one loaded, else `ambiguous`. A destructive operation skips the one being worked
 * with: with several loaded, it acts only on a named structure.
 */
export function oneStructure(
  ctx: ApiContext,
  name: string | undefined,
  options: { readonly destructive?: boolean } = {},
): StructureRef | ApiFailure {
  const list = loadedStructures(ctx);
  if (!name && list.length > 1 && !options.destructive) {
    const working = workingStructure(ctx.plugin, list);
    if (working) return working.structure;
  }
  return pickNamed(list, name, structureRule(list));
}

/** The named structure, or every loaded one. */
export function someStructures(ctx: ApiContext, name: string | undefined): readonly StructureRef[] | ApiFailure {
  const list = loadedStructures(ctx);
  if (!list.length) return ApiResult.fail('not-found', 'Nothing is loaded.', 'Load a structure first.');
  if (!name) return list;
  const s = pickNamed(list, name, structureRule(list));
  return isFailure(s) ? s : [s];
}

// --- components -------------------------------------------------------------------------------------------------

/** The component's name for people. */
export function componentName(c: StructureComponentRef): string {
  const params = c.cell.params?.values as StructureComponentParams | undefined;
  return params?.label || c.cell.obj?.label || c.key || c.cell.transform.ref;
}

export interface FoundComponent {
  readonly structure: StructureRef;
  readonly component: StructureComponentRef;
}

/** One component among these structures' (by state ref or key, else by name), by the one rule. */
export function findComponent(
  ctx: ApiContext,
  structures: readonly StructureRef[],
  name: string | undefined,
): FoundComponent | ApiFailure {
  const all = structures.flatMap((structure) => structure.components.map((component) => ({ structure, component })));
  return pickNamed(all, name, {
    what: 'component',
    ids: (x) =>
      x.component.key ? [x.component.cell.transform.ref, x.component.key] : [x.component.cell.transform.ref],
    names: (x) => [componentName(x.component)],
    describe: (x) =>
      `${componentName(x.component)} (${x.component.cell.transform.ref}${structures.length > 1 ? `, in ${labelOf(ctx, x.structure)}` : ''})`,
    none: {
      message: 'There are no components.',
      hint: structures.length === 1 ? `${labelOf(ctx, structures[0])} has none.` : undefined,
    },
    listed: 'Components',
    several: 'components are there',
    more: structures.length > 1 ? ' or name the structure' : '',
  });
}

/** One component by the one rule, in the named structure or in all of them (a component ref is unique among them). */
export function inStructures(ctx: ApiContext, name: string, structure: string): FoundComponent | ApiFailure {
  const structures = someStructures(ctx, structure || undefined);
  return isFailure(structures) ? structures : findComponent(ctx, structures, name);
}

// --- hierarchy lists (volumes, particle lists) ------------------------------------------------------------------

type Listed = { cell: { obj?: { label?: string }; transform: { ref: string } } };

/** Names for what is in a hierarchy list (volumes, particle lists): the label, told apart by the ref when it repeats. */
export function namesOf(list: readonly Listed[]): ReadonlyMap<string, string> {
  const counts = new Map<string, number>();
  for (const x of list) counts.set(x.cell.obj?.label ?? '', (counts.get(x.cell.obj?.label ?? '') ?? 0) + 1);
  const names = new Map<string, string>();
  for (const x of list) {
    const label = x.cell.obj?.label || x.cell.transform.ref;
    names.set(
      x.cell.transform.ref,
      (counts.get(x.cell.obj?.label ?? '') ?? 0) > 1 ? `${label} · ${x.cell.transform.ref}` : label,
    );
  }
  return names;
}

/** One of a hierarchy list by the one rule; `what` is "volume" or "particle list", `empty` the hint when none is loaded. */
export function pickOne<T extends Listed>(
  list: readonly T[],
  name: string | undefined,
  what: string,
  empty: string,
): T | ApiFailure {
  const names = namesOf(list);
  const shown = (x: T) => names.get(x.cell.transform.ref) ?? x.cell.transform.ref;
  return pickNamed(list, name, {
    what,
    ids: (x) => [x.cell.transform.ref],
    names: (x) => [shown(x), x.cell.obj?.label],
    describe: (x) => (shown(x) === x.cell.transform.ref ? shown(x) : `${shown(x)} (${x.cell.transform.ref})`),
    none: { message: `No ${what}s are loaded.`, hint: empty },
  });
}

// --- the state tree ---------------------------------------------------------------------------------------------

/** The highest ancestor that holds only this node (the model, trajectory, data or format it came from), below the root. */
export function entryRoot(ctx: ApiContext, ref: string): string {
  const tree = ctx.plugin.state.data.tree;
  let current = ref;
  for (;;) {
    const parent = tree.transforms.get(current)?.parent;
    if (!parent || parent === StateTransform.RootRef || tree.children.get(parent).size !== 1) return current;
    current = parent;
  }
}
