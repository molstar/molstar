/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type State, type StateObjectCell, StateTransform } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { setSubtreeVisibility } from '@molstar/plugin/behavior/static/state';
import { type ApiAreaContext, ApiArea, type ApiContext, ApiResult } from '../../api/area.js';
import { actionName } from '../../api/load.js';
import { mergeParams } from '../../api/params.js';
import { checkValues } from '../../api/run.js';
import { chainBase, createStateTreeSource, foldChain, type StateCellNode, type StateTreeSource } from './source.js';

// `state` (operations.md §2, the `state-tree` feature): the plugin's state tree for advanced use. `get` reads it; `remove`,
// `visibility`, `update` and `action` change one cell by its state ref, through `plugin.state.data` as the classic tree
// does. Cells the other areas manage (a structure, a representation) are the same cells: removing a structure here is
// what `structure.remove` does to its structure. Each is a history step (the run's checkpoint restores it).

const RefParam = () => PD.Text('', { description: 'A state cell by ref: see state.get' });

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

function cellOf(ctx: ApiContext, ref: string): StateObjectCell | ReturnType<typeof ApiResult.fail> {
  if (!ref) return ApiResult.fail('invalid', 'Which cell?', 'Give its state ref: state.get lists them.');
  const cell = ctx.plugin.state.data.cells.get(ref);
  if (!cell) return ApiResult.fail('not-found', `There is no state cell '${ref}'.`, 'state.get lists the cells.');
  return cell;
}

const labelOf = (cell: StateObjectCell) => cell.obj?.label ?? cell.transform.transformer.definition.display.name;

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
      ref: PD.Text('', { description: 'One cell, with its params and applicable actions; empty: the whole tree' }),
    },
    read: true,
    run: async ({ ref }, ctx): Promise<ApiResult<StateCellInfo[] | StateCellDetails>> => {
      const source = this.tree;
      if (!ref) {
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
      const cell = cellOf(ctx, ref);
      if (ApiResult.isFailure(cell)) return cell;
      const details = source.details(ref);
      if (!details) return ApiResult.fail('not-found', `There is no state cell '${ref}'.`);
      const { children: _, ...node } = details.node;
      const own = details.updates.find((u) => u.ref === node.ref);
      const value: StateCellDetails = {
        ...info(node as StateCellNode),
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
        `${node.label} (${node.ref}) is ${node.transform}${node.visible === false ? ', hidden' : ''}. Params: ${brief(value.params)}.${actions}`,
        {},
        { value },
      );
    },
  });

  readonly visibility = this.operation({
    id: 'state.visibility',
    title: 'Show or hide a state cell',
    params: { ref: RefParam(), visible: PD.Boolean(true) },
    async run({ ref, visible }, ctx) {
      const cell = cellOf(ctx, ref);
      if (ApiResult.isFailure(cell)) return cell;
      if (ref === StateTransform.RootRef) return ApiResult.fail('blocked', 'The root has nothing to show or hide.');
      const name = labelOf(cell);
      if (!!cell.state.isHidden === !visible) {
        return ApiResult.already(`${name} is already ${visible ? 'shown' : 'hidden'}.`);
      }
      // The cell and what hangs off it, as the classic tree's eye.
      setSubtreeVisibility(ctx.plugin.state.data, ref, !visible);
      return ApiResult.ok(`${visible ? 'Showed' : 'Hid'} ${name}.`, { nodes: [ref] });
    },
  });

  readonly update = this.operation({
    id: 'state.update',
    title: 'Update a state cell',
    surface: { kind: 'agent-only', reason: "The State tree tab's forms." },
    params: {
      ref: RefParam(),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The transform's params to change; the others stay as they are" },
      ),
    },
    async run({ ref, params }, ctx) {
      const cell = cellOf(ctx, ref);
      if (ApiResult.isFailure(cell)) return cell;
      const name = labelOf(cell);
      if (!cell.params || ref === StateTransform.RootRef)
        return ApiResult.fail('blocked', `${name} has no params to update.`);
      // The patch against the transformer's params: unknown names, wrong types, a select's options, a number's range.
      const checked = checkValues(cell.params.definition, params);
      if (ApiResult.isFailure(checked))
        return ApiResult.fail(checked.code, `${name}: ${checked.message}`, checked.hint);
      const current = cell.transform.params as Record<string, any>;
      const next = mergeParams(current, params);
      if (PD.areEqual(cell.params.definition, current, next)) return ApiResult.already(`${name} is already so.`);
      const state = ctx.plugin.state.data;
      const failedBefore = new Set(erroredIn(state, ref).map((c) => c.transform.ref));
      await state.build().to(ref).update(next).commit();
      // The cell or something built on it failed with the new params: the run's checkpoint restores the old ones.
      const failed = erroredIn(state, ref).find((c) => !failedBefore.has(c.transform.ref));
      if (failed) return ApiResult.fail('failed', `${name}: ${failed.errorText ?? 'the update failed'}.`);
      return ApiResult.ok(`Updated ${name}.`, { nodes: [ref] });
    },
  });

  readonly action = this.operation({
    id: 'state.action',
    title: 'Apply an action to a state cell',
    surface: { kind: 'agent-only', reason: "The State tree tab's forms." },
    params: {
      ref: RefParam(),
      action: PD.Text('', { description: 'An action applicable to the cell, by name: add-representation, …' }),
      params: PD.Value<Record<string, unknown>>(
        {},
        { description: "The action's params; missing ones are its defaults" },
      ),
    },
    async run({ ref, action, params }, ctx) {
      const state = ctx.plugin.state.data;
      const cell = cellOf(ctx, ref);
      if (ApiResult.isFailure(cell)) return cell;
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
    title: 'Remove a state cell',
    params: { ref: RefParam() },
    approval: 'destructive',
    async run({ ref }, ctx) {
      const cell = cellOf(ctx, ref);
      if (ApiResult.isFailure(cell)) return cell;
      const state = ctx.plugin.state.data;
      if (ref === StateTransform.RootRef)
        return ApiResult.fail('blocked', 'The root cannot be removed.', 'Clear the session instead.');
      const name = labelOf(cell);
      if (cell.state.isLocked) return ApiResult.fail('blocked', `${name} is locked.`);
      const count = subtreeSize(state, ref);
      await state.build().delete(removalRoot(state, ref)).commit();
      return ApiResult.ok(`Removed ${name}${count > 1 ? ` and ${count - 1} below it` : ''}.`, { nodes: [ref] });
    },
  });
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
