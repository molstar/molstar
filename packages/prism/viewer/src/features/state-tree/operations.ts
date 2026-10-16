/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type State, type StateObjectCell, StateTransform, StateTransformer } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { setSubtreeVisibility } from '@molstar/plugin/behavior/static/state';
import { type ApiAreaContext, ApiArea, type ApiContext, type ApiFailure, ApiResult } from '../../api/area.js';
import { type NodeKind, nodeKind, nodesFor, selectNodes, selectOneNode } from '../../api/nodes.js';
import { actionName } from '../../api/load.js';
import { mergeParams, NodeParam } from '../../api/params.js';
import { checkValues } from '../../api/run.js';
import { chainBase, createStateTreeSource, foldChain, type StateCellNode, type StateTreeSource } from './source.js';

// `state` (operations.md §2, the `state-tree` feature): the plugin's state tree for advanced use. `get` reads it and
// `select` lists what a node selector matches (api/nodes.ts); `visibility`, `update` and `remove` change every cell the
// selector matches, `action` one, and `apply` adds a transform under each, through `plugin.state.data` as the classic
// tree does. Cells the other areas manage (a structure, a representation) are the same cells: removing a structure
// here is what `structure.remove` does to its structure. Each change is a history step (the run's checkpoint restores
// it).

const CellsParam = () => NodeParam('State cells by selector (structure > component, #ref) or a ref: see state.get');
const CellParam = () => NodeParam('One state cell by selector or ref: see state.get');

/** What `state.select` gives for a cell. */
export interface StateNodeInfo {
  readonly ref: string;
  readonly label: string;
  readonly kind?: NodeKind;
  /** The transformer's id. */
  readonly transformer: string;
  /** The parent's ref. */
  readonly parent: string;
  readonly visible: boolean;
  readonly status: 'ok' | 'pending' | 'error';
}

/** What `state.get` gives for a cell. */
export interface StateCellInfo {
  readonly ref: string;
  readonly label: string;
  readonly description?: string;
  /** The transformer's id. */
  readonly transformer: string;
  /** The transformer's display name. */
  readonly transform: string;
  readonly decorators?: readonly string[];
  readonly visible?: boolean;
  readonly status: 'ok' | 'pending' | 'error';
  readonly error?: string;
  readonly children?: readonly StateCellInfo[];
}

/** One cell in full: its params as values, the updates it takes and the actions applicable to it. */
export interface StateCellDetails extends Omit<StateCellInfo, 'children'> {
  readonly params: Record<string, unknown>;
  /** The cell's own transform and the decorators folded into it, each with its params. */
  readonly updates: readonly {
    readonly ref: string;
    readonly name: string;
    readonly params: Record<string, unknown>;
  }[];
  /** The cell actions apply to. */
  readonly target: string;
  readonly actions: readonly {
    readonly name: string;
    readonly title: string;
    readonly description?: string;
    /** The names of the params it takes (defaults fill the others). */
    readonly params: readonly string[];
  }[];
}

/** Lines at most, so a big session does not fill a model's context; the value has all of it. */
const OutlineLines = 200;

const info = (n: StateCellNode): StateCellInfo => ({
  ref: n.ref,
  label: n.label,
  ...(n.description ? { description: n.description } : {}),
  transformer: n.transformer,
  transform: n.transform,
  ...(n.decorators ? { decorators: n.decorators } : {}),
  ...(n.visible !== undefined ? { visible: n.visible } : {}),
  status: n.status ?? 'ok',
  ...(n.error ? { error: n.error } : {}),
  ...(n.children ? { children: n.children.map(info) } : {}),
});

function outline(nodes: readonly StateCellNode[], lines: string[], depth = 0): void {
  for (const n of nodes) {
    if (lines.length > OutlineLines) return;
    const flags = [
      n.visible === false ? 'hidden' : '',
      n.status === 'error' ? `error: ${n.error ?? 'failed'}` : n.status === 'pending' ? 'pending' : '',
    ].filter(Boolean);
    lines.push(
      `${'  '.repeat(depth)}${n.ref}: ${n.label}${n.description ? ` (${n.description})` : ''} [${n.transform}]${flags.length ? ` ${flags.join(', ')}` : ''}`,
    );
    if (n.children) outline(n.children, lines, depth + 1);
  }
}

/** A value as JSON for a result line, shortened (a param may hold a long list). */
function brief(value: unknown, max = 400): string {
  try {
    const text = JSON.stringify(value) ?? '';
    return text.length > max ? `${text.slice(0, max)}…` : text;
  } catch {
    return '(not shown)';
  }
}

function cellOf(ctx: ApiContext, node: string): StateObjectCell | ApiFailure {
  if (!node) return ApiResult.fail('invalid', 'Which cell?', 'Give a node selector or a ref: state.get lists them.');
  return selectOneNode(ctx.plugin.state.data, node);
}

function cellsOf(ctx: ApiContext, node: string): StateObjectCell[] | ApiFailure {
  return nodesFor(ctx.plugin.state.data, node, 'state cell', (cell) => cell);
}

/** "Cartoon", "3 cells". */
const cellsName = (cells: readonly StateObjectCell[]) =>
  cells.length === 1 ? labelOf(cells[0]) : `${cells.length} cells`;

/** A transformer by id, or by its name without the namespace (`structure-representation-3d`). */
function transformerOf(name: string): StateTransformer | undefined {
  if (StateTransformer.has(name)) return StateTransformer.get(name);
  const found = StateTransformer.getAll().filter((t) => t.id.endsWith(`.${name}`));
  return found.length === 1 ? found[0] : undefined;
}

const labelOf = (cell: StateObjectCell) => cell.obj?.label ?? cell.transform.transformer.definition.display.name;

const words = (key: string) => key.replace(/([a-z0-9])([A-Z])/g, '$1 $2').toLowerCase();

/** "Cartoon of 1CBS": the cell and the structure it belongs to (the cell's own name for a structure or above). */
function nameWithStructure(ctx: ApiContext, cell: StateObjectCell): string {
  const name = labelOf(cell);
  const state = ctx.plugin.state.data;
  const known = new Map(ctx.sources.structures.value.map((s) => [s.ref, s.label]));
  let ref = cell.transform.ref;
  while (ref !== StateTransform.RootRef) {
    const structure = known.get(ref);
    if (structure) return ref === cell.transform.ref ? structure : `${name} of ${structure}`;
    const parent = state.transforms.get(ref)?.parent;
    if (!parent || parent === ref) break;
    ref = parent;
  }
  return name;
}

/** The params that differ, in words ("color theme", "size theme"); a long list is counted. */
function changedParams(params: PD.Params, before: Record<string, any>, after: Record<string, any>): string {
  const keys = Object.keys(params).filter((k) => !PD.isParamEqual(params[k], before[k], after[k]));
  const names = keys.map((k) => (params[k].label ?? words(k)).toLowerCase());
  return names.length > 3 ? `${names.length} params` : names.join(', ');
}

/** The ancestors that are ghosts with this cell as their only child go with it (the classic tree's Remove). */
function removalRoot(state: State, ref: string): string {
  const tree = state.tree;
  let curr = tree.transforms.get(ref);
  while (curr.parent !== curr.ref) {
    if (tree.children.get(curr.parent).size > 1) return curr.ref;
    const parent = tree.transforms.get(curr.parent);
    if (!parent.state.isGhost) return curr.ref;
    curr = parent;
  }
  return curr.ref;
}

export class StateApi extends ApiArea {
  readonly id = 'state';

  /** The tree for views and view models: followed only while someone subscribes (`source.ts`). */
  readonly tree: StateTreeSource;

  constructor(ctx: ApiAreaContext) {
    super(ctx);
    this.tree = createStateTreeSource(ctx.plugin);
  }

  dispose(): void {
    this.tree.dispose();
  }

  readonly get = this.operation({
    id: 'state.get',
    title: 'Read the state tree',
    params: {
      node: NodeParam('One cell (selector or ref), with its params and applicable actions; empty: the whole tree'),
    },
    read: true,
    run: async ({ node }, ctx): Promise<ApiResult<StateCellInfo[] | StateCellDetails>> => {
      const source = this.tree;
      if (!node) {
        const roots = source.roots.value;
        const lines: string[] = [];
        outline(roots, lines);
        const more = lines.length > OutlineLines ? `\n… (${OutlineLines} lines of it; the value has all)` : '';
        return ApiResult.ok(
          lines.length ? `${lines.slice(0, OutlineLines).join('\n')}${more}` : 'The state is empty.',
          {},
          { value: roots.map(info) },
        );
      }
      const cell = cellOf(ctx, node);
      if (ApiResult.isFailure(cell)) return cell;
      const ref = cell.transform.ref;
      const details = source.details(ref);
      if (!details) return ApiResult.fail('not-found', `There is no state cell '${ref}'.`);
      const { children: _, ...cellNode } = details.node;
      const own = details.updates.find((u) => u.ref === cellNode.ref);
      const value: StateCellDetails = {
        ...info(cellNode as StateCellNode),
        params: (own?.values ?? {}) as Record<string, unknown>,
        updates: details.updates.map((u) => ({ ref: u.ref, name: u.name, params: u.values })),
        target: details.target,
        actions: details.actions.map((a) => ({
          name: a.name,
          title: a.title,
          ...(a.description ? { description: a.description } : {}),
          params: Object.keys(a.params),
        })),
      };
      const actions = value.actions.length ? ` Actions: ${value.actions.map((a) => a.name).join(', ')}.` : '';
      return ApiResult.ok(
        `${cellNode.label} (${cellNode.ref}) is ${cellNode.transform}${cellNode.visible === false ? ', hidden' : ''}. Params: ${brief(value.params)}.${actions}`,
        {},
        { value },
      );
    },
  });

  readonly visibility = this.operation({
    id: 'state.visibility',
    title: 'Show or hide state cells',
    params: { node: CellsParam(), visible: PD.Boolean(true) },
    async run({ node, visible }, ctx) {
      const cells = cellsOf(ctx, node);
      if (ApiResult.isFailure(cells)) return cells;
      const shown = cells.filter((c) => c.transform.ref !== StateTransform.RootRef);
      if (!shown.length) return ApiResult.fail('blocked', 'The root has nothing to show or hide.');
      const todo = shown.filter((c) => !!c.state.isHidden === visible);
      if (!todo.length) {
        return ApiResult.already(
          `${cellsName(shown)} ${shown.length > 1 ? 'are' : 'is'} already ${visible ? 'shown' : 'hidden'}.`,
        );
      }
      // Each cell and what hangs off it, as the classic tree's eye.
      for (const c of todo) setSubtreeVisibility(ctx.plugin.state.data, c.transform.ref, !visible);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${cellsName(todo)}.`, {
        nodes: todo.map((c) => c.transform.ref),
      });
    },
  });

  readonly update = this.operation({
    id: 'state.update',
    title: 'Update state cells',
    surface: { kind: 'agent-only', reason: "The State tree tab's forms." },
    params: {
      node: CellsParam(),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The transform's params to change; the others stay as they are" },
      ),
    },
    async run({ node, params }, ctx) {
      const cells = cellsOf(ctx, node);
      if (ApiResult.isFailure(cells)) return cells;
      const state = ctx.plugin.state.data;
      // Each is checked before any changes, so a bad patch changes nothing.
      const plans: { cell: StateObjectCell; current: Record<string, any>; next: Record<string, any> }[] = [];
      for (const cell of cells) {
        const name = labelOf(cell);
        if (!cell.params || cell.transform.ref === StateTransform.RootRef)
          return ApiResult.fail('blocked', `${name} has no params to update.`);
        // Where data comes from (the root's children: a download, a file) is the loaders' and their approvals' (api/area.ts).
        if (cell.transform.parent === StateTransform.RootRef) {
          return ApiResult.fail(
            'blocked',
            `${name} is where data comes from.`,
            'Load other data with load.entry or load.url.',
          );
        }
        // The patch against the transformer's params: unknown names, wrong types, a select's options, a number's range.
        const checked = checkValues(cell.params.definition, params);
        if (ApiResult.isFailure(checked))
          return ApiResult.fail(checked.code, `${name}: ${checked.message}`, checked.hint);
        const current = cell.transform.params as Record<string, any>;
        const next = mergeParams(current, params);
        if (!PD.areEqual(cell.params.definition, current, next)) plans.push({ cell, current, next });
      }
      if (!plans.length) return ApiResult.already(`${cellsName(cells)} ${cells.length > 1 ? 'are' : 'is'} already so.`);
      const failedBefore = new Set(cells.flatMap((c) => erroredIn(state, c.transform.ref)).map((c) => c.transform.ref));
      const b = state.build();
      for (const { cell, next } of plans) b.to(cell.transform.ref).update(next);
      await b.commit();
      // A cell or something built on it failed with the new params: the run's checkpoint restores the old ones.
      for (const { cell } of plans) {
        const failed = erroredIn(state, cell.transform.ref).find((c) => !failedBefore.has(c.transform.ref));
        if (failed) return ApiResult.fail('failed', `${labelOf(cell)}: ${failed.errorText ?? 'the update failed'}.`);
      }
      const message =
        plans.length === 1
          ? (() => {
              const { cell, current, next } = plans[0];
              const what = changedParams(cell.params!.definition, current, next);
              return `Updated ${nameWithStructure(ctx, cell)}${what ? `: ${what}` : ''}.`;
            })()
          : `Updated ${plans.length} cells.`;
      return ApiResult.ok(message, { nodes: plans.map((p) => p.cell.transform.ref) });
    },
  });

  readonly action = this.operation({
    id: 'state.action',
    title: 'Apply an action to a state cell',
    // Actions download (Download Structure, Download Density on the root): network (api/area.ts).
    approval: 'network',
    surface: { kind: 'agent-only', reason: "The State tree tab's forms." },
    params: {
      node: CellParam(),
      action: PD.Text('', { description: 'An action applicable to the cell, by name: add-representation, …' }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The action's params; missing ones are its defaults" },
      ),
    },
    async run({ node, action, params }, ctx) {
      const state = ctx.plugin.state.data;
      const cell = cellOf(ctx, node);
      if (ApiResult.isFailure(cell)) return cell;
      const ref = cell.transform.ref;
      // Actions go to the end of a decorator chain, which the cell's children hang off (the classic tree's too).
      const chain = foldChain(state, chainBase(state, ref));
      const target = chain[chain.length - 1] ?? cell;
      const actions = state.actions.fromCell(target, ctx.plugin);
      const found = actions.find((a) => actionName(a) === action || a.id === action);
      if (!found) {
        return ApiResult.fail(
          'not-found',
          action
            ? `${labelOf(cell)} has no action '${action}'.`
            : `Which action? ${labelOf(cell)} takes ${actions.length || 'none'}.`,
          actions.length
            ? `Actions: ${actions.map(actionName).join(', ')}.`
            : 'No action applies to it; state.get with its ref says so.',
        );
      }
      const name = found.definition.display.name;
      const defs = (target.obj && found.definition.params?.(target.obj, ctx.plugin)) as PD.Params | undefined;
      const checked = checkValues(defs ?? {}, params);
      if (ApiResult.isFailure(checked))
        return ApiResult.fail(checked.code, `${name}: ${checked.message}`, checked.hint);
      try {
        await ctx.plugin.runTask(state.applyAction(found, checked.values as never, target.transform.ref));
      } catch (e) {
        return ApiResult.fail('failed', `${name}: ${e instanceof Error ? e.message : String(e)}`);
      }
      return ApiResult.ok(`${name}: done.`, { nodes: [target.transform.ref] });
    },
  });

  readonly remove = this.operation({
    id: 'state.remove',
    title: 'Remove state cells',
    params: { node: CellsParam() },
    approval: 'destructive',
    async run({ node }, ctx) {
      const cells = cellsOf(ctx, node);
      if (ApiResult.isFailure(cells)) return cells;
      const state = ctx.plugin.state.data;
      if (cells.some((c) => c.transform.ref === StateTransform.RootRef))
        return ApiResult.fail('blocked', 'The root cannot be removed.', 'Clear the session instead.');
      const locked = cells.find((c) => c.state.isLocked);
      if (locked) return ApiResult.fail('blocked', `${labelOf(locked)} is locked.`);
      // A cell below another that goes, goes with it.
      const refs = new Set(cells.map((c) => c.transform.ref));
      const tops = cells.filter((c) => !hasAncestorIn(state, c.transform.ref, refs));
      const count = tops.reduce((n, c) => n + subtreeSize(state, c.transform.ref), 0);
      const b = state.build();
      for (const c of tops) b.delete(removalRoot(state, c.transform.ref));
      await b.commit();
      const below = count - tops.length;
      return ApiResult.ok(`Removed ${cellsName(tops)}${below > 0 ? ` and ${below} below` : ''}.`, {
        nodes: tops.map((c) => c.transform.ref),
      });
    },
  });

  // operations.md §2: what a node selector matches, as a list (the tree's reading of a selector; nothing changes).
  readonly select = this.operation({
    id: 'state.select',
    title: 'Find state cells',
    params: { node: NodeParam('A node selector: structure > component[type=ligand] representation, :hidden, #ref') },
    read: true,
    run: async ({ node }, ctx): Promise<ApiResult<StateNodeInfo[]>> => {
      const state = ctx.plugin.state.data;
      const cells = selectNodes(state, node || '*');
      if (ApiResult.isFailure(cells)) return cells;
      const value = cells.map((c): StateNodeInfo => {
        const kind = nodeKind(state, c);
        return {
          ref: c.transform.ref,
          label: labelOf(c),
          ...(kind ? { kind } : {}),
          transformer: c.transform.transformer.id,
          parent: c.transform.parent,
          visible: !c.state.isHidden,
          status: c.status === 'error' ? 'error' : c.status === 'ok' ? 'ok' : 'pending',
        };
      });
      if (!value.length) return ApiResult.ok(`Nothing matches '${node}'.`, {}, { value });
      const lines = value
        .slice(0, OutlineLines)
        .map((v) => `${v.ref}: ${v.label}${v.kind ? ` (${v.kind})` : ''}${v.visible ? '' : ', hidden'}`);
      const more = value.length > OutlineLines ? `\n… (${OutlineLines} of ${value.length}; the value has all)` : '';
      return ApiResult.ok(`${value.length} match:\n${lines.join('\n')}${more}`, {}, { value });
    },
  });

  // operations.md §2: a transform applied under every cell a selector matches (Mol*'s `build().to(ref).apply(...)`);
  // the new cells' refs are the value. One history step: undo removes them.
  readonly apply = this.operation({
    id: 'state.apply',
    title: 'Apply a transform',
    surface: { kind: 'agent-only', reason: 'Scripts and the agent; the panels add what they show through the areas.' },
    params: {
      node: CellsParam(),
      transformer: PD.Text('', {
        description: 'A transformer by id or short name: structure-representation-3d, structure-component, …',
      }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The transform's params; missing ones are its defaults" },
      ),
      tags: PD.Value<string[]>([], { description: 'Tags for the new cells, to select them by later (.tag)' }),
    },
    run: async ({ node, transformer, params, tags }, ctx): Promise<ApiResult<{ refs: string[] }>> => {
      const cells = cellsOf(ctx, node);
      if (ApiResult.isFailure(cells)) return cells;
      const t = transformerOf(transformer);
      if (!t) {
        const known = cells[0]?.obj ? StateTransformer.fromType(cells[0].obj.type).map((x) => x.id) : [];
        return ApiResult.fail(
          transformer ? 'not-found' : 'invalid',
          transformer ? `There is no transformer '${transformer}' (or several end so).` : 'Which transformer?',
          known.length
            ? `For ${labelOf(cells[0])}: ${known.slice(0, 20).join(', ')}${known.length > 20 ? ', …' : ''}.`
            : undefined,
        );
      }
      const name = t.definition.display.name;
      const state = ctx.plugin.state.data;
      const b = state.build();
      const refs: string[] = [];
      for (const cell of cells) {
        // Where data comes from is the loaders' (api/area.ts): no downloads under the root here.
        if (cell.transform.ref === StateTransform.RootRef) {
          return ApiResult.fail('blocked', `${name} cannot go on the root.`, 'Load data with load.entry or load.url.');
        }
        if (!cell.obj || !t.definition.from.some((f) => f.type === cell.obj!.type)) {
          return ApiResult.fail(
            'invalid',
            `${name} does not apply to ${labelOf(cell)} (${cell.obj?.type.name ?? 'not ready'}).`,
            cell.obj
              ? `For it: ${StateTransformer.fromType(cell.obj.type)
                  .map((x) => x.id)
                  .slice(0, 20)
                  .join(', ')}.`
              : undefined,
          );
        }
        if (t.definition.isApplicable && !t.definition.isApplicable(cell.obj, ctx.plugin)) {
          return ApiResult.fail('invalid', `${name} does not apply to ${labelOf(cell)}.`);
        }
        const defs = StateTransformer.getParamDefinition(t, cell.obj, ctx.plugin) as PD.Params;
        const checked = checkValues(defs, params);
        if (ApiResult.isFailure(checked))
          return ApiResult.fail(checked.code, `${name}: ${checked.message}`, checked.hint);
        const added = b.to(cell.transform.ref).apply(t, checked.values as never, tags.length ? { tags } : undefined);
        refs.push(added.ref);
      }
      await b.commit();
      const failed = refs.map((r) => state.cells.get(r)).find((c) => !c || c.status === 'error');
      if (failed) return ApiResult.fail('failed', `${name}: ${failed.errorText ?? 'it failed'}.`);
      return ApiResult.ok(`Applied ${name} to ${cellsName(cells)}.`, { nodes: refs }, { value: { refs } });
    },
  });
}

function hasAncestorIn(state: State, ref: string, refs: ReadonlySet<string>): boolean {
  let parent = state.tree.transforms.get(ref)?.parent;
  while (parent && parent !== ref) {
    if (refs.has(parent)) return true;
    const next = state.tree.transforms.get(parent)?.parent;
    if (!next || next === parent) return false;
    ref = parent;
    parent = next;
  }
  return false;
}

/** The cells of a subtree (itself included) that are in error. */
function erroredIn(state: State, ref: string): StateObjectCell[] {
  const out: StateObjectCell[] = [];
  const walk = (r: string) => {
    const cell = state.cells.get(r);
    if (cell?.status === 'error') out.push(cell);
    for (const c of state.tree.children.get(r).toArray()) walk(c);
  };
  walk(ref);
  return out;
}

function subtreeSize(state: State, ref: string): number {
  let n = 0;
  const walk = (r: string) => {
    n++;
    for (const c of state.tree.children.get(r).toArray()) walk(c);
  };
  walk(ref);
  return n;
}
