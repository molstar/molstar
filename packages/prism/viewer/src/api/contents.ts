/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { StructureComponent } from '@molstar/plugin/state/transforms/structure/selection';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import { StructureElement, StructureProperties } from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
import type { PluginContext } from '@molstar/plugin/context';
import type { StructureRef } from '@molstar/plugin/state/manager/structure/hierarchy-state';
import { all } from '@molstar/plugin/state/queries/structure/basic';
import { nonStandardPolymer } from '@molstar/plugin/state/queries/structure/residue';
import {
  branchedPlusConnected,
  coarse,
  ion,
  ligandPlusConnected,
  lipid,
  nucleic,
  polymer,
  protein,
  water,
} from '@molstar/plugin/state/queries/structure/type';
import { type ApiAreaContext, type ApiContext, ApiArea, ApiResult } from './area.js';
import { EveryStructureParam, StructureParam } from './params.js';
import { isFailure, labelOf, type StructureId, structureId, loadedStructures, oneStructure } from './scope.js';

// `contents` (api.md §2; README decision 18): hiding part of a structure at the entity and chain level, so what
// is hidden is left out of every component (and so of every representation) of the structure, not just made
// transparent. Each component is made to hold "what it held, except what is hidden"; its own params are kept in the
// viewer's state (`registerCustomState`, so snapshots, sessions and the history's checkpoints carry them) and put back
// when nothing is hidden any more. Components made later (a preset, `component.add`) get the same (`reapplyContents`).

/** The viewer's state for this, by structure ref: what is hidden, and each changed component's own params. */
interface ContentsState {
  readonly version: 1;
  readonly hidden: Readonly<
    Record<string, { readonly entities: readonly string[]; readonly chains: readonly string[] }>
  >;
  readonly originals: Readonly<Record<string, unknown>>;
}

const StateKey = 'prism-viewer.contents';
const Empty: ContentsState = { version: 1, hidden: {}, originals: {} };

/** The expression of a static component type (`structure-component.ts`'s queries). */
const StaticQueries: Readonly<Record<string, Expression>> = {
  all: all.expression,
  polymer: polymer.expression,
  protein: protein.expression,
  nucleic: nucleic.expression,
  water: water.expression,
  ion: ion.expression,
  lipid: lipid.expression,
  branched: branchedPlusConnected.expression,
  ligand: ligandPlusConnected.expression,
  'non-standard': nonStandardPolymer.expression,
  coarse: coarse.expression,
};

/** One per viewer: the area holds it, the scene's other operations reapply it. */
const states = new WeakMap<object, { value: ContentsState }>();

/** What is hidden of a structure (the scene's contents rows show it); empty when nothing is. */
export function hiddenContents(
  plugin: object,
  structureRef: string,
): { readonly entities: readonly string[]; readonly chains: readonly string[] } {
  return states.get(plugin)?.value.hidden[structureRef] ?? NothingHidden;
}

/** What `hiddenContents` says when nothing is hidden: one value, so a reader can tell "as before" by identity. */
const NothingHidden: { readonly entities: readonly string[]; readonly chains: readonly string[] } = Object.freeze({
  entities: Object.freeze([]),
  chains: Object.freeze([]),
});

/** A component left empty by hidden contents: its ref, the label it had, and its parent (where it stands). */
export interface EmptiedComponent {
  readonly ref: string;
  readonly label: string;
  readonly parent: string;
}

/**
 * The components of each structure in `structures` that hold nothing while its contents are hidden (their cells are
 * empty, so the hierarchy leaves them out): the scene keeps their rows. With each, the label it had. Read once for all
 * the structures: the saved params walked once (per structure, a sync over S structures was S × the saved params;
 * AGENTS.md, bulk operations). A component is under every structure among its ancestors, in the saved params' order.
 */
export function emptiedComponentsOf(
  plugin: PluginContext,
  structures: ReadonlySet<string>,
): ReadonlyMap<string, readonly EmptiedComponent[]> {
  const out = new Map<string, EmptiedComponent[]>();
  const originals = states.get(plugin)?.value.originals;
  if (!originals) return out;
  const data = plugin.state.data;
  for (const ref in originals) {
    const own = originals[ref] as ComponentParams;
    const cell = data.cells.get(ref);
    if (!cell || cell.transform.transformer !== StructureComponent) continue;
    const held = cell.obj?.data as { elementCount?: number } | undefined;
    if (held?.elementCount) continue;
    let e: EmptiedComponent | undefined;
    // Its ancestors (itself first, as `under` reads them).
    let current: string | undefined = ref;
    while (current) {
      if (structures.has(current)) {
        e ??= {
          ref,
          label: own.label || (cell.transform.params as { label?: string } | undefined)?.label || 'Component',
          parent: cell.transform.parent,
        };
        let list = out.get(current);
        if (!list) out.set(current, (list = []));
        list.push(e);
      }
      const parent: string | undefined = data.tree.transforms.get(current)?.parent;
      if (!parent || parent === current) break;
      current = parent;
    }
  }
  return out;
}

/** Whether a cell is in a subtree (by the tree's parents). */
function under(tree: { transforms: { get(ref: string): { parent: string } | undefined } }, ref: string, root: string) {
  let current: string | undefined = ref;
  while (current) {
    if (current === root) return true;
    const parent: string | undefined = tree.transforms.get(current)?.parent;
    if (!parent || parent === current) return false;
    current = parent;
  }
  return false;
}

function stateOf(ctx: ApiAreaContext): { value: ContentsState } {
  let s = states.get(ctx.plugin);
  if (!s) {
    s = { value: Empty };
    states.set(ctx.plugin, s);
  }
  return s;
}

export class ContentsApi extends ApiArea {
  readonly id = 'contents';
  readonly #unregister: () => void;

  constructor(ctx: ApiAreaContext) {
    super(ctx);
    const state = stateOf(ctx);
    this.#unregister = ctx.plugin.state.registerCustomState(StateKey, {
      get: () => (Object.keys(state.value.hidden).length ? state.value : undefined),
      set: (value) => {
        const v = value as ContentsState | undefined;
        state.value = v?.version === 1 ? v : Empty;
      },
    });
  }

  readonly visibility = this.operation({
    id: 'contents.visibility',
    title: 'Show or hide an entity or a chain',
    params: {
      entity: PD.Text('', { description: "An entity by id ('1', '2', …); with chain, that chain of it only" }),
      chain: PD.Text('', {
        description:
          "A chain by author id ('A', 'B', …): every entity's (a ligand's and waters' too), unless entity is given",
      }),
      visible: PD.Boolean(true),
      structure: StructureParam(),
    },
    async run({ entity, chain, visible, structure }, ctx) {
      if (!entity && !chain)
        return ApiResult.fail('invalid', 'Give an entity or a chain.', "entity: '1', chain: 'A', or both.");
      const s = oneStructure(ctx, structure || undefined);
      if (isFailure(s)) return s;
      const ref = s.cell.transform.ref;
      const data = s.cell.obj?.data;
      if (!data) return ApiResult.fail('not-found', `${labelOf(ctx, s)} has no data.`);
      const known = contentsOf(s);
      if (entity && !known.entities.has(entity)) {
        return ApiResult.fail(
          'not-found',
          `${labelOf(ctx, s)} has no entity '${entity}'.`,
          `Entities: ${[...known.entities].join(', ')}.`,
        );
      }
      if (entity && chain && !known.pairs.has(`${entity}:${chain}`)) {
        return ApiResult.fail(
          'not-found',
          `Entity ${entity} of ${labelOf(ctx, s)} has no chain '${chain}'.`,
          `Its chains: ${[...known.pairs]
            .filter((p) => p.startsWith(`${entity}:`))
            .map((p) => p.slice(entity.length + 1))
            .join(', ')}.`,
        );
      }
      if (chain && !known.chains.has(chain)) {
        return ApiResult.fail(
          'not-found',
          `${labelOf(ctx, s)} has no chain '${chain}'.`,
          `Chains: ${[...known.chains].join(', ')}.`,
        );
      }
      const state = stateOf(ctx);
      const current = state.value.hidden[ref] ?? { entities: [], chains: [] };
      // A chain of an entity is kept as `entity:chain` (only its atoms); a chain alone, every entity's in it.
      const key = chain ? 'chains' : 'entities';
      const id = entity && chain ? `${entity}:${chain}` : entity || chain;
      const isHidden = current[key].includes(id);
      const what =
        entity && chain ? `chain ${chain} of entity ${entity}` : entity ? `entity ${entity}` : `chain ${chain}`;
      if (isHidden === !visible) {
        return ApiResult.already(
          `${what[0].toUpperCase()}${what.slice(1)} of ${labelOf(ctx, s)} is already ${visible ? 'shown' : 'hidden'}.`,
        );
      }
      const next = { ...current, [key]: visible ? current[key].filter((x) => x !== id) : [...current[key], id] };
      const hidden = { ...state.value.hidden };
      if (next.entities.length || next.chains.length) hidden[ref] = next;
      else delete hidden[ref];
      state.value = { ...state.value, hidden };
      await reapplyContents(ctx, [s]);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${what} of ${labelOf(ctx, s)}.`, {
        nodes: s.components.map((c) => c.cell.transform.ref),
      });
    },
  });

  readonly get = this.operation({
    id: 'contents.get',
    title: 'List the contents',
    params: { structure: EveryStructureParam() },
    read: true,
    async run({ structure }, ctx) {
      const list = structure ? [oneStructure(ctx, structure)] : loadedStructures(ctx);
      const lines: string[] = [];
      const value: {
        structure: StructureId;
        entities: string[];
        chains: string[];
        hidden: { entities: readonly string[]; chains: readonly string[] };
      }[] = [];
      for (const s of list) {
        if (isFailure(s)) return s;
        const c = contentsOf(s);
        const hidden = stateOf(ctx).value.hidden[s.cell.transform.ref] ?? { entities: [], chains: [] };
        value.push({ structure: structureId(ctx, s), entities: [...c.entities], chains: [...c.chains], hidden });
        const h = [...hidden.entities.map((e) => `entity ${e}`), ...hidden.chains.map((x) => `chain ${x}`)];
        lines.push(
          `${labelOf(ctx, s)}: entities ${[...c.entities].join(', ')}; chains ${[...c.chains].join(', ')}${h.length ? `; hidden: ${h.join(', ')}` : ''}`,
        );
      }
      return ApiResult.ok(lines.join('\n') || 'Nothing is loaded.', {}, { value });
    },
  });

  dispose(): void {
    this.#unregister();
  }
}

/** A structure's entities and chains (author ids), and which chains each entity has (`entity:chain`). */
function contentsOf(s: StructureRef): { entities: Set<string>; chains: Set<string>; pairs: Set<string> } {
  const entities = new Set<string>();
  const chains = new Set<string>();
  const pairs = new Set<string>();
  const data = s.cell.obj?.data;
  if (!data) return { entities, chains, pairs };
  const loc = StructureElement.Location.create(data);
  for (const unit of data.units) {
    loc.unit = unit;
    loc.element = unit.elements[0];
    const entity = StructureProperties.chain.label_entity_id(loc);
    const chain = StructureProperties.chain.auth_asym_id(loc);
    entities.add(entity);
    chains.add(chain);
    pairs.add(`${entity}:${chain}`);
  }
  return { entities, chains, pairs };
}

/**
 * Makes the structures' components leave out what is hidden in them (and hold what they held when nothing is). Called by
 * `contents.visibility`, and by the operations that make components (a preset, `component.add`), inside their step.
 */
export async function reapplyContents(ctx: ApiContext, structures: readonly StructureRef[]): Promise<void> {
  const state = stateOf(ctx);
  const data = ctx.plugin.state.data;
  const originals = { ...state.value.originals };
  const build = data.build();
  let changed = false;
  for (const s of structures) {
    const structureRef = s.cell.transform.ref;
    const hidden = state.value.hidden[structureRef];
    // A component that holds nothing once the hidden part is left out is a null cell, not in the hierarchy: the
    // state tree has it, so it comes back when shown again.
    const components = new Set(s.components.map((c) => c.cell.transform.ref));
    for (const ref of Object.keys(originals)) if (isUnder(ctx, ref, structureRef)) components.add(ref);
    for (const ref of components) {
      const cell = data.cells.get(ref);
      // Components only: the focus's target and surroundings are in the hierarchy's list too.
      if (!cell || cell.transform.transformer !== StructureComponent) continue;
      if (!hidden) {
        if (ref in originals) {
          build.to(ref).update(originals[ref] as never);
          delete originals[ref];
          changed = true;
        }
        continue;
      }
      const own = (originals[ref] ?? cell.transform.params) as ComponentParams;
      const kept = expressionOf(own);
      if (!kept) continue;
      originals[ref] = own;
      build.to(ref).update({
        ...own,
        label: own.label || cell.obj?.label || '',
        type: { name: 'expression', params: MS.struct.modifier.exceptBy({ 0: kept, by: hiddenExpression(hidden) }) },
      } as never);
      changed = true;
    }
  }
  // Components whose cells are gone keep no originals.
  for (const ref of Object.keys(originals)) if (!data.cells.has(ref)) delete originals[ref];
  state.value = { ...state.value, originals };
  if (changed) await build.commit();
}

/** A component's own params while its structure hides contents (its cell holds them minus what is hidden). */
export function ownParams(ctx: ApiAreaContext, ref: string): ComponentParams | undefined {
  return stateOf(ctx).value.originals[ref] as ComponentParams | undefined;
}

/**
 * Sets a component's own params (a rename, a modification): while its structure hides contents they are kept as its
 * own and the cell gets them minus what is hidden, so showing the contents again keeps the change.
 */
export async function setOwnParams(
  ctx: ApiContext,
  s: StructureRef,
  ref: string,
  params: ComponentParams,
): Promise<void> {
  const state = stateOf(ctx);
  if (!(ref in state.value.originals)) {
    await ctx.plugin.state.data
      .build()
      .to(ref)
      .update(params as never)
      .commit();
    return;
  }
  state.value = { ...state.value, originals: { ...state.value.originals, [ref]: params } };
  await reapplyContents(ctx, [s]);
}

/** Whether a cell is in a structure's subtree. */
function isUnder(ctx: ApiContext, ref: string, root: string): boolean {
  return under(ctx.plugin.state.data.tree, ref, root);
}

export interface ComponentParams {
  type:
    | { name: 'static'; params: string }
    | { name: 'expression'; params: Expression }
    | { name: 'bundle'; params: StructureElement.Bundle }
    | { name: 'script'; params: Script };
  label?: string;
  nullIfEmpty?: boolean;
}

export function expressionOf(p: ComponentParams): Expression | undefined {
  switch (p.type.name) {
    case 'static':
      return StaticQueries[p.type.params];
    case 'expression':
      return p.type.params;
    case 'bundle':
      return StructureElement.Bundle.toExpression(p.type.params);
    case 'script':
      return Script.toExpression(p.type.params);
  }
}

function hiddenExpression(hidden: {
  readonly entities: readonly string[];
  readonly chains: readonly string[];
}): Expression {
  const tests: Expression[] = [];
  if (hidden.entities.length) {
    tests.push(MS.core.set.has([MS.set(...hidden.entities), MS.ammp('label_entity_id')]));
  }
  const whole = hidden.chains.filter((c) => !c.includes(':'));
  if (whole.length) tests.push(MS.core.set.has([MS.set(...whole), MS.ammp('auth_asym_id')]));
  // A chain of an entity: its atoms only (a ligand or waters with the same author chain stay).
  for (const pair of hidden.chains) {
    const at = pair.indexOf(':');
    if (at < 0) continue;
    tests.push(
      MS.core.logic.and([
        MS.core.rel.eq([MS.ammp('label_entity_id'), pair.slice(0, at)]),
        MS.core.rel.eq([MS.ammp('auth_asym_id'), pair.slice(at + 1)]),
      ]),
    );
  }
  return MS.struct.generator.atomGroups({ 'chain-test': tests.length === 1 ? tests[0] : MS.core.logic.or(tests) });
}
