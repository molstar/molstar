/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { MolScriptBuilder as MS } from '@molstar/query-language/language/builder';
import type { Expression } from '@molstar/query-language/language/expression';
import { StructureElement, StructureProperties } from '@molstar/model/model/structure';
import { Script } from '@molstar/model/script/script';
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
import { isFailure, labelOf, loadedStructures, oneStructure } from './scope.js';

// `contents` (operations.md §2; README decision 18): hiding part of a structure at the entity and chain level, so what
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
      entity: PD.Text('', { description: "An entity by id ('1', '2', …)" }),
      chain: PD.Text('', { description: "A chain by author id ('A', 'B', …)" }),
      visible: PD.Boolean(true),
      structure: PD.Text('', { description: 'A structure by label or state ref; empty: the one in view' }),
    },
    async run({ entity, chain, visible, structure }, ctx) {
      if (!entity === !chain)
        return ApiResult.fail('invalid', 'Give an entity or a chain.', "entity: '1', or chain: 'A'.");
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
      if (chain && !known.chains.has(chain)) {
        return ApiResult.fail(
          'not-found',
          `${labelOf(ctx, s)} has no chain '${chain}'.`,
          `Chains: ${[...known.chains].join(', ')}.`,
        );
      }
      const state = stateOf(ctx);
      const current = state.value.hidden[ref] ?? { entities: [], chains: [] };
      const key = entity ? 'entities' : 'chains';
      const id = entity || chain;
      const isHidden = current[key].includes(id);
      const what = entity ? `entity ${entity}` : `chain ${chain}`;
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
    params: { structure: PD.Text('', { description: 'A structure by label or state ref; empty: every one' }) },
    read: true,
    async run({ structure }, ctx) {
      const list = structure ? [oneStructure(ctx, structure)] : loadedStructures(ctx);
      const lines: string[] = [];
      const value: {
        structure: string;
        entities: string[];
        chains: string[];
        hidden: { entities: readonly string[]; chains: readonly string[] };
      }[] = [];
      for (const s of list) {
        if (isFailure(s)) return s;
        const c = contentsOf(s);
        const hidden = stateOf(ctx).value.hidden[s.cell.transform.ref] ?? { entities: [], chains: [] };
        value.push({ structure: labelOf(ctx, s), entities: [...c.entities], chains: [...c.chains], hidden });
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

/** A structure's entities and chains (author ids). */
function contentsOf(s: StructureRef): { entities: Set<string>; chains: Set<string> } {
  const entities = new Set<string>();
  const chains = new Set<string>();
  const data = s.cell.obj?.data;
  if (!data) return { entities, chains };
  const loc = StructureElement.Location.create(data);
  for (const unit of data.units) {
    loc.unit = unit;
    loc.element = unit.elements[0];
    entities.add(StructureProperties.chain.label_entity_id(loc));
    chains.add(StructureProperties.chain.auth_asym_id(loc));
  }
  return { entities, chains };
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
      if (!cell) continue;
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

/** Whether a cell is in a structure's subtree. */
function isUnder(ctx: ApiContext, ref: string, root: string): boolean {
  const tree = ctx.plugin.state.data.tree;
  let current: string | undefined = ref;
  while (current) {
    if (current === root) return true;
    const parent: string | undefined = tree.transforms.get(current)?.parent;
    if (!parent || parent === current) return false;
    current = parent;
  }
  return false;
}

interface ComponentParams {
  type:
    | { name: 'static'; params: string }
    | { name: 'expression'; params: Expression }
    | { name: 'bundle'; params: StructureElement.Bundle }
    | { name: 'script'; params: Script };
  label?: string;
  nullIfEmpty?: boolean;
}

function expressionOf(p: ComponentParams): Expression | undefined {
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
  if (hidden.chains.length) {
    tests.push(MS.core.set.has([MS.set(...hidden.chains), MS.ammp('auth_asym_id')]));
  }
  return MS.struct.generator.atomGroups({ 'chain-test': tests.length === 1 ? tests[0] : MS.core.logic.or(tests) });
}
