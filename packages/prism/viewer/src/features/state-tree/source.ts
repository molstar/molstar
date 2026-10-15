/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type State, type StateObjectCell, StateObject, StateTransform } from '@molstar/core/state';
import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject as SO } from '@molstar/plugin/state/objects';
import { type IconName, observableState, type ReadonlyState } from '@molstar/prism-ui/core';
import { Observable, Subject, Subscription } from 'rxjs';
import { actionName } from '../../api/load.js';
import type { StateTreeNode } from '../../model/state-tree.js';
import { sameItems } from '../../sources/state.js';

// The state tree as the viewer reads it (viewer-spec.md §3.1, the `state-tree` feature): the cells of
// `plugin.state.data` as nodes (decorators folded into the cell they decorate, ghost cells passed through, as the
// classic tree shows them), and the details of one cell: its updates and applicable actions. DOM-free, and with
// `operations.ts` the only code of the feature that reads the plugin (`_test/boundary.test.ts`).
//
// The tree follows `state.data.events` while someone looks at it: events only mark cells, one flush per short delay
// rebuilds the marked cells and their ancestors, and every other node keeps its identity (a keyed list re-renders
// what is a different value, so nothing else is touched).

/** A node of the tree: the view's node, and which transformer made it. */
export interface StateCellNode extends Omit<StateTreeNode, 'children'> {
  /** The transformer's id ('ms-plugin.download'); `transform` is its display name. */
  readonly transformer: string;
  readonly children?: readonly StateCellNode[];
}

/** One of a cell's transforms that can be updated: its own, then the decorators folded into it. */
export interface CellUpdate {
  /** The cell the update goes to. */
  readonly ref: string;
  /** The transformer's display name. */
  readonly name: string;
  readonly params: PD.Params;
  readonly values: Record<string, unknown>;
}

/** A state action applicable to a cell. */
export interface CellAction {
  /** The display name as a slug (`add-representation`): the plugin's action ids are made per session. */
  readonly name: string;
  readonly title: string;
  readonly description?: string;
  /** What running it does ("Add", "Download"). */
  readonly trigger?: string;
  readonly params: PD.Params;
}

export interface CellDetails {
  readonly node: StateCellNode;
  /** The cell actions apply to: the end of the decorator chain, which the cell's children hang off. */
  readonly target: string;
  readonly updates: readonly CellUpdate[];
  readonly actions: readonly CellAction[];
}

export interface StateTreeSource {
  /** The tree. Followed only while subscribed (a hidden tab holds no subscription to the plugin). */
  readonly roots: ReadonlyState<readonly StateCellNode[]>;
  /** Cells whose params (or object) changed or went, after each flush: the forms over them reload. */
  readonly changed: Observable<ReadonlySet<string>>;
  /** The node of a cell that has one. */
  node(ref: string): StateCellNode | undefined;
  /** A cell by its ref, or any cell of a decorator chain by the chain's. */
  details(ref: string): CellDetails | undefined;
  /** Rebuilds what events marked now (tests; the delay is for coalescing). */
  flush(): void;
  dispose(): void;
}

/** A cell and the decorators folded into it (A → A transforms that follow it, each the only child of the one before). */
export function foldChain(state: State, ref: string): StateObjectCell[] {
  const first = state.cells.get(ref);
  if (!first) return [];
  const chain = [first];
  let end = first;
  while (end.status === 'ok') {
    const kids = state.tree.children.get(end.transform.ref);
    if (kids.size !== 1) break;
    const next = state.cells.get(kids.first()!);
    if (!next || next.status !== 'ok' || !next.transform.transformer.definition.isDecorator) break;
    chain.push(next);
    end = next;
  }
  return chain;
}

/** The cell a decorator belongs to, else the cell itself. */
export function chainBase(state: State, ref: string): string {
  let cell = state.cells.get(ref);
  while (cell?.transform.transformer.definition.isDecorator && cell.transform.parent !== cell.transform.ref) {
    const parent = state.cells.get(cell.transform.parent);
    if (!parent) break;
    cell = parent;
  }
  return cell?.transform.ref ?? ref;
}

const displayName = (cell: StateObjectCell) => cell.transform.transformer.definition.display.name;

const Visuals = new Set(['cartoon', 'ball-and-stick', 'spacefill']);

function iconOf(chain: readonly StateObjectCell[]): IconName {
  const base = chain[0];
  const end = chain[chain.length - 1];
  const obj = end.obj;
  if (obj && SO.isRepresentation3D(obj)) {
    const type = (end.transform.params as { type?: { name?: string } } | undefined)?.type?.name ?? '';
    if (Visuals.has(type)) return type as IconName;
    if (type.includes('surface')) return 'surface';
    return 'shapes';
  }
  const name = displayName(base);
  const typeName = obj?.type.name ?? '';
  if (name === 'Download') return 'download';
  if (name === 'Structure Component') return 'layers';
  if (typeName === 'Structure') return 'structure';
  if (typeName === 'Trajectory') return 'trajectory';
  if (typeName === 'Model') return 'box';
  if (typeName.includes('Volume')) return 'volume';
  if (typeName === 'CIF File') return 'file-code';
  if (obj?.type.typeClass === 'Data') return 'file';
  if (obj?.type.typeClass === 'Group' || obj?.type.typeClass === 'Root') return 'layers';
  return 'shapes';
}

/** The node of a cell and its chain; `children` are what hangs off the chain's end. */
function makeNode(chain: readonly StateObjectCell[], children: readonly StateCellNode[] | undefined): StateCellNode {
  const base = chain[0];
  const end = chain[chain.length - 1];
  const obj = end.obj;
  const status = end.status === 'error' ? 'error' : end.status === 'ok' ? 'ok' : 'pending';
  const hasVisual = (!!obj && SO.isRepresentation3D(obj)) || !!children?.some((c) => c.visible !== undefined);
  const description = obj?.description;
  return {
    ref: base.transform.ref,
    label: obj && end.status === 'ok' ? obj.label : displayName(end),
    ...(description ? { description } : {}),
    icon: iconOf(chain),
    transform: displayName(base),
    transformer: base.transform.transformer.id,
    ...(chain.length > 1 ? { decorators: chain.slice(1).map(displayName) } : {}),
    ...(children?.length ? { children } : {}),
    ...(hasVisual ? { visible: !base.state.isHidden } : {}),
    status,
    ...(status === 'error' && end.errorText ? { error: end.errorText } : {}),
  };
}

function sameNode(a: StateCellNode, b: StateCellNode): boolean {
  return (
    a.ref === b.ref &&
    a.label === b.label &&
    a.description === b.description &&
    a.icon === b.icon &&
    a.transform === b.transform &&
    a.transformer === b.transformer &&
    a.visible === b.visible &&
    a.status === b.status &&
    a.error === b.error &&
    sameItems(a.decorators ?? [], b.decorators ?? []) &&
    sameItems(a.children ?? [], b.children ?? [])
  );
}

export interface StateTreeSourceOptions {
  /** How long events are coalesced (ms); default 40. */
  readonly delay?: number;
}

export function createStateTreeSource(plugin: PluginContext, options: StateTreeSourceOptions = {}): StateTreeSource {
  const state = plugin.state.data;
  const delay = options.delay ?? 40;
  const roots = observableState<readonly StateCellNode[]>([], { compare: sameItems });
  const changed = new Subject<ReadonlySet<string>>();
  // What each cell shows (by the ref of the cell that starts a chain): none (null object), one node, or, for a ghost
  // cell, its children's.
  const cache = new Map<string, readonly StateCellNode[]>();
  let invalid: ReadonlySet<string> | 'all' = 'all';

  const build = (ref: string): readonly StateCellNode[] => {
    const cell = state.cells.get(ref);
    if (!cell || cell.obj === StateObject.Null || !state.transforms.has(ref)) return [];
    if (invalid !== 'all' && !invalid.has(ref)) {
      const cached = cache.get(ref);
      if (cached) return cached;
    }
    const old = cache.get(ref);
    const chain = foldChain(state, ref);
    const children = childrenOf(chain[chain.length - 1].transform.ref);
    let made: readonly StateCellNode[];
    if (cell.state.isGhost && cell.status === 'ok') made = children;
    else {
      const node = makeNode(chain, children);
      made = old?.length === 1 && old[0].ref === ref && sameNode(old[0], node) ? old : [node];
    }
    if (old && sameItems(old, made)) made = old;
    cache.set(ref, made);
    return made;
  };

  const childrenOf = (ref: string): StateCellNode[] => {
    const out: StateCellNode[] = [];
    for (const child of state.tree.children.get(ref).toArray()) out.push(...build(child));
    return out;
  };

  // --- events -> marks -> one flush ---------------------------------------------------------------------------

  const dirty = new Set<string>();
  const paramsDirty = new Set<string>();
  let all = true;
  let removed = false;
  let timer: ReturnType<typeof setTimeout> | undefined;

  const schedule = () => {
    if (timer === undefined) timer = setTimeout(flush, delay);
  };
  const mark = (ref: string | undefined, params = false) => {
    if (ref === undefined) return;
    dirty.add(ref);
    if (params) paramsDirty.add(ref);
    schedule();
  };

  function flush() {
    if (timer !== undefined) clearTimeout(timer);
    timer = undefined;
    if (!all && !dirty.size && !paramsDirty.size) return;
    if (all) invalid = 'all';
    else {
      // A changed cell and every ancestor: the ancestors' children lists hold it.
      const set = new Set<string>();
      for (const ref of dirty) {
        let r: string | undefined = ref;
        while (r !== undefined && !set.has(r)) {
          set.add(r);
          const t = state.transforms.get(r);
          r = t && t.parent !== r ? t.parent : undefined;
        }
      }
      invalid = set;
    }
    const params = new Set(paramsDirty);
    const prune = all || removed;
    dirty.clear();
    paramsDirty.clear();
    all = false;
    removed = false;
    roots.set(childrenOf(StateTransform.RootRef));
    invalid = new Set();
    if (prune) for (const ref of [...cache.keys()]) if (!state.cells.has(ref)) cache.delete(ref);
    if (params.size) changed.next(params);
  }

  const events = state.events;
  const listen = (): Subscription => {
    const subs = new Subscription();
    subs.add(
      events.cell.created.subscribe((e) => {
        mark(e.ref);
        mark(e.cell.transform.parent);
      }),
    );
    subs.add(
      events.cell.removed.subscribe((e) => {
        removed = true;
        mark(e.parent);
        paramsDirty.add(e.ref);
      }),
    );
    subs.add(events.cell.stateUpdated.subscribe((e) => mark(e.ref)));
    subs.add(events.object.created.subscribe((e) => mark(e.ref, true)));
    subs.add(events.object.updated.subscribe((e) => mark(e.ref, true)));
    subs.add(events.object.removed.subscribe((e) => mark(e.ref, true)));
    // The end of an update that restructured the tree without a cell event (a restored snapshot).
    subs.add(events.changed.subscribe(() => schedule()));
    return subs;
  };

  // --- activity: the plugin is followed while someone subscribes ------------------------------------------------

  let watchers = 0;
  let listening: Subscription | undefined;
  const activate = () => {
    if (watchers++ > 0) return;
    // Whatever happened while nobody looked: the whole tree again (nodes that are the same keep their identity).
    all = true;
    flush();
    listening = listen();
  };
  const deactivate = () => {
    if (--watchers > 0) return;
    listening?.unsubscribe();
    listening = undefined;
    if (timer !== undefined) clearTimeout(timer);
    timer = undefined;
    dirty.clear();
    paramsDirty.clear();
    all = true;
  };
  const watched = <T>(source: Observable<T>) =>
    new Observable<T>((subscriber) => {
      activate();
      const sub = source.subscribe(subscriber);
      return () => {
        sub.unsubscribe();
        deactivate();
      };
    });
  const current = () => {
    if (!listening) {
      all = true;
      flush();
    }
    return roots.value;
  };

  const rootsState = watched(roots.asReadonly) as ReadonlyState<readonly StateCellNode[]>;
  Object.defineProperty(rootsState, 'value', { get: current });

  const nodeOf = (ref: string) => {
    current();
    return cache.get(ref)?.find((n) => n.ref === ref);
  };

  return {
    roots: rootsState,
    changed: watched(changed.asObservable()),
    node: nodeOf,
    details(ref) {
      const base = chainBase(state, ref);
      const chain = foldChain(state, base);
      if (!chain.length) return undefined;
      const end = chain[chain.length - 1];
      const node = nodeOf(base) ?? makeNode(chain, undefined);
      const updates = chain
        .filter((c) => !!c.params)
        .map(
          (c): CellUpdate => ({
            ref: c.transform.ref,
            name: displayName(c),
            params: c.params!.definition,
            values: c.params!.values,
          }),
        );
      const actions = [...state.actions.fromCell(end, plugin)]
        .sort((a, b) => a.definition.display.name.localeCompare(b.definition.display.name))
        .map((a): CellAction => {
          let params: PD.Params = {};
          try {
            params = (end.obj && a.definition.params?.(end.obj, plugin)) || {};
          } catch {
            // An action whose params cannot be made for this object is not offered with a form.
          }
          const { description, triggerLabel } = a.definition.display;
          return {
            name: actionName(a),
            title: a.definition.display.name,
            ...(description ? { description } : {}),
            ...(triggerLabel ? { trigger: triggerLabel } : {}),
            params,
          };
        });
      return { node, target: end.transform.ref, updates, actions };
    },
    flush,
    dispose() {
      listening?.unsubscribe();
      listening = undefined;
      if (timer !== undefined) clearTimeout(timer);
      timer = undefined;
      changed.complete();
    },
  };
}
