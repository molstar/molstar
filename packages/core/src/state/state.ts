/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Taylor Hoffmann <taylor@hoffmann.io>
 */

import { StateObject, type StateObjectCell, StateObjectSelector } from './object.js';
import { StateTree } from './tree.js';
import { StateTransform } from './transform.js';
import { StateTransformer } from './transformer.js';
import { type RuntimeContext, Task } from '@molstar/core/task';
import { StateSelection } from './state/selection.js';
import { RxEventHelper } from '@molstar/core/util/rx-event-helper';
import { StateBuilder } from './state/builder.js';
import { findDeletes, findEditedDeletes, findEditedRoots, findUpdateRoots } from './state/diff.js';
import type { StateAction } from './action.js';
import { StateActionManager } from './action/manager.js';
import type { TransientTree } from './tree/transient.js';
import { LogEntry } from '@molstar/core/util/log-entry';
import { now, formatTimespan } from '@molstar/core/util/now';
import { ParamDefinition } from '@molstar/core/util/param-definition';
import { StateTreeSpine } from './tree/spine.js';
import { AsyncQueue } from '@molstar/core/util/async-queue';
import { arraySetAdd, arraySetRemove } from '@molstar/core/util/array';
import { UniqueArray } from '@molstar/core/data/generic/unique-array';
import { assignIfUndefined } from '@molstar/core/util/object';
import { SortedArray } from '@molstar/core/data/int/sorted-array';

export { State, StateTreeCycleError };

/**
 * Thrown when a cycle is detected in the state-tree dependency graph
 * (the cross-edges defined by `transform.dependsOn` and effective
 * dependencies derived from params). `cycle` holds the closed path,
 * e.g. `['A', 'B', 'A']`.
 */
class StateTreeCycleError extends Error {
  readonly cycle: StateTransform.Ref[];
  constructor(cycle: StateTransform.Ref[]) {
    super(`Cyclic state-tree dependency detected: ${cycle.join(' -> ')}`);
    this.name = 'StateTreeCycleError';
    this.cycle = cycle;
  }
}

class State {
  private _tree: TransientTree;

  protected errorFree = true;
  /**
   * the cells are in line with the tree (README invariant 2): false after an update that threw or a fast track with an
   * uncounted write, until a whole-tree diff completes (reconcile.md §2)
   */
  protected cellsInSync = true;

  private ev = RxEventHelper.create();

  readonly globalContext: unknown = void 0;
  readonly events = {
    cell: {
      stateUpdated: this.ev<State.ObjectEvent & { cell: StateObjectCell }>(),
      created: this.ev<State.ObjectEvent & { cell: StateObjectCell }>(),
      removed: this.ev<State.ObjectEvent & { parent: StateTransform.Ref }>(),
    },
    object: {
      updated: this.ev<
        State.ObjectEvent & { action: 'in-place' | 'recreate'; obj: StateObject; oldObj?: StateObject; oldData?: any }
      >(),
      created: this.ev<State.ObjectEvent & { obj: StateObject }>(),
      removed: this.ev<State.ObjectEvent & { obj?: StateObject }>(),
    },
    log: this.ev<LogEntry>(),
    changed: this.ev<{ state: State; inTransaction: boolean }>(),
    historyUpdated: this.ev<{ state: State }>(),
  };

  readonly behaviors = {
    currentObject: this.ev.behavior<State.ObjectEvent>({ state: this, ref: StateTransform.RootRef }),
    isUpdating: this.ev.behavior<boolean>(false),
  };

  readonly actions = new StateActionManager();

  readonly runTask: <T>(task: Task<T>) => Promise<T>;

  get tree(): StateTree {
    return this._tree;
  }
  get transforms() {
    return (this._tree as StateTree).transforms;
  }
  get current() {
    return this.behaviors.currentObject.value.ref;
  }
  get root() {
    return this.cells.get((this._tree as StateTree).root.ref)!;
  }

  build() {
    const builder = new StateBuilder.Root(this.tree, this);
    this.builderFlagMarks.set(builder.editInfo, this.flagWrites);
    return builder;
  }

  /**
   * Flag writes by `updateCellState`: ref → flag → the number of its latest write and the value written. A builder's
   * tree is a snapshot of the tree it was made on, which keeps its identity across flag writes, so its commit carries
   * the flags written later over (tree.md §4), and only those: a flag the builder set itself stays.
   */
  private flagWrites = 0;
  private flagLog: FlagLog = new Map();
  private builderFlagMarks = new WeakMap<StateBuilder.EditInfo, number>();
  /**
   * The flag writes in order (write number, ref), so a commit visits only the refs written after its builder's mark,
   * not every ref ever flagged (perf-audit-2 03/P17). Compacted to each ref's last write when it grows.
   */
  private flagWriteOrder: { n: number[]; refs: StateTransform.Ref[] } = { n: [], refs: [] };

  private logFlagWrite(n: number, ref: StateTransform.Ref) {
    const order = this.flagWriteOrder;
    order.n.push(n);
    order.refs.push(ref);
    if (order.n.length > 2 * this.flagLog.size + 64) {
      // a ref is written after a mark iff its last write is: keep the last write of each ref still logged
      const last: [number, StateTransform.Ref][] = [];
      for (const [r, writes] of this.flagLog) {
        let m = 0;
        for (const k of StateFlags) if (writes[k] && writes[k]!.n > m) m = writes[k]!.n;
        last.push([m, r]);
      }
      last.sort((a, b) => a[0] - b[0]);
      this.flagWriteOrder = { n: last.map((e) => e[0]), refs: last.map((e) => e[1]) };
    }
  }

  private carryFlagWrites(editInfo: StateBuilder.EditInfo, tree: TransientTree) {
    const mark = this.builderFlagMarks.get(editInfo) ?? 0;
    if (mark === this.flagWrites) return;
    const order = this.flagWriteOrder;
    const visited = new Set<StateTransform.Ref>();
    for (let i = SortedArray.findPredecessorIndex(order.n as any as SortedArray, mark + 1); i < order.n.length; i++) {
      const ref = order.refs[i];
      if (visited.has(ref)) continue;
      visited.add(ref);
      const writes = this.flagLog.get(ref);
      if (!writes || !tree.transforms.has(ref)) continue;
      const state = tree.transforms.get(ref).state;
      let update: Partial<StateTransform.State> | undefined;
      for (const k of StateFlags) {
        const w = writes[k];
        if (!w || w.n <= mark || !!state[k] === !!w.value) continue;
        if (!update) update = {};
        update[k] = w.value;
      }
      if (update) tree.assignState(ref, update);
    }
  }

  readonly cells: State.Cells = new Map();
  private spine = new StateTreeSpine.Impl(this.cells);

  private refResolvers: [name: string, resolver: (state: State, ref: StateTransform.Ref) => StateObject | undefined][] =
    [];

  registerRefResolver(name: string, resolver: (state: State, ref: StateTransform.Ref) => StateObject | undefined) {
    this.refResolvers.push([name, resolver]);
  }

  removeRefResolver(name: string) {
    this.refResolvers = this.refResolvers.filter((r) => r[0] !== name);
  }

  tryGetCellData = <T extends StateObject>(ref: StateTransform.Ref) => {
    let ret = this.cells.get(ref)?.obj?.data;
    if (ret === undefined) {
      for (const [, resolver] of this.refResolvers) {
        ret = resolver(this, ref);
        if (ret !== undefined) break;
      }
      if (ret === undefined) throw new Error(`Cell '${ref}' data undefined.`);
    }
    return ret as T extends StateObject<infer D> ? D : never;
  };

  private historyCapacity = 5;
  private history: [StateTree, string][] = [];

  private addHistory(tree: StateTree, label?: string) {
    if (this.historyCapacity === 0) return;

    this.history.unshift([tree, label || 'Update']);
    if (this.history.length > this.historyCapacity) this.history.pop();

    this.events.historyUpdated.next({ state: this });
  }

  private clearHistory() {
    if (this.history.length === 0) return;
    this.history = [];
    this.events.historyUpdated.next({ state: this });
  }

  get latestUndoLabel() {
    return this.history.length > 0 ? this.history[0][1] : void 0;
  }

  get canUndo() {
    return this.history.length > 0;
  }

  private undoingHistory = false;

  undo() {
    return Task.create('Undo', async (ctx) => {
      const e = this.history.shift();
      if (!e) return;
      this.events.historyUpdated.next({ state: this });
      this.undoingHistory = true;
      try {
        await this.updateTree(e[0], { canUndo: false }).runInContext(ctx);
      } finally {
        this.undoingHistory = false;
      }
    });
  }

  getSnapshot(): State.Snapshot {
    return { tree: StateTree.toJSON(this._tree) };
  }

  setSnapshot(snapshot: State.Snapshot, options?: Partial<State.UpdateOptions>) {
    const tree = StateTree.fromJSON(snapshot.tree);
    return this.updateTree(tree, options);
  }

  setCurrent(ref: StateTransform.Ref) {
    this.behaviors.currentObject.next({ state: this, ref });
  }

  updateCellState(
    ref: StateTransform.Ref,
    stateOrProvider: ((old: StateTransform.State) => Partial<StateTransform.State>) | Partial<StateTransform.State>,
  ) {
    const cell = this.cells.get(ref);
    if (!cell) return;

    const update = typeof stateOrProvider === 'function' ? stateOrProvider(cell.state) : stateOrProvider;

    const written = update && StateFlags.filter((k) => k in update && !!update[k] !== !!cell.state[k]);
    if (StateTransform.assignState(cell.state, update)) {
      cell.transform = this._tree.assignState(cell.transform.ref, update);
      const n = ++this.flagWrites;
      let writes = this.flagLog.get(ref);
      if (!writes) this.flagLog.set(ref, (writes = {}));
      for (const k of written!) writes[k] = { n, value: update![k] };
      this.logFlagWrite(n, ref);
      this.events.cell.stateUpdated.next({ state: this, ref, cell });
    }
  }

  dispose() {
    // Dispose every still-live cell so transformer dispose callbacks
    // (e.g. WebGL/GL buffer cleanup) actually run. Without this,
    // calling dispose() on a State that still has cells leaks any
    // resources held by transformer dispose callbacks because they
    // would only fire on per-cell deletion (see updateNode/findDeletes).
    const refs: StateTransform.Ref[] = [];
    StateTree.doPostOrder(this._tree, this._tree.root, { refs }, (n, _, s) => {
      s.refs.push(n.ref);
    });
    for (let i = refs.length - 1; i >= 0; i--) {
      const cell = (this.cells as Map<StateTransform.Ref, StateObjectCell>).get(refs[i]);
      if (!cell) continue;
      try {
        dispose(cell.transform, cell.obj, cell.transform.params, cell.cache, this.globalContext);
      } catch (e) {
        console.warn('Error in transformer dispose during State.dispose', e);
      }
    }

    this.ev.dispose();
    this.actions.dispose();
  }

  /**
   * Select Cells using the provided selector.
   * @example state.query(StateSelection.Generators.byRef('test').ancestorOfType(type))
   * @example state.query('test')
   */
  select<C extends StateObjectCell>(selector: StateSelection.Selector<C>) {
    return StateSelection.select(selector, this);
  }

  /**
   * Select Cells by building a query generated on the fly.
   * @example state.select(q => q.byRef('test').subtree())
   */
  selectQ<C extends StateObjectCell>(selector: (q: typeof StateSelection.Generators) => StateSelection.Selector<C>) {
    if (typeof selector === 'string') return StateSelection.select(selector, this);
    return StateSelection.select(selector(StateSelection.Generators), this);
  }

  /**
   * Creates a Task that applies the specified StateAction (i.e. must use run* on the result)
   * If no ref is specified, apply to root.
   */
  applyAction<A extends StateAction>(
    action: A,
    params: StateAction.Params<A>,
    ref: StateTransform.Ref = StateTransform.RootRef,
  ): Task<void> {
    return Task.create('Apply Action', (ctx) => {
      const cell = this.cells.get(ref);
      if (!cell) throw new Error(`'${ref}' does not exist.`);
      if (cell.status !== 'ok') throw new Error(`Action cannot be applied to a cell with status '${cell.status}'`);

      return runTask(action.definition.run({ ref, cell, a: cell.obj!, params, state: this }, this.globalContext), ctx);
    });
  }

  private inTransaction = false;
  private inTransactionError = false;

  /** Apply series of updates to the state. If any of them fail, revert to the original state. */
  transaction(
    edits: (ctx: RuntimeContext) => Promise<void> | void,
    options?: { canUndo?: string | boolean; rethrowErrors?: boolean },
  ) {
    return Task.create('State Transaction', async (ctx) => {
      const isNested = this.inTransaction;

      // if (!isNested) this.changedInTransaction = false;

      const snapshot = this._tree.asImmutable();
      let restored = false;
      try {
        if (!isNested) this.behaviors.isUpdating.next(true);

        this.inTransaction = true;
        this.inTransactionError = false;
        await edits(ctx);

        if (this.inTransactionError) {
          restored = true;
          await this.updateTree(snapshot).runInContext(ctx);
        }
      } catch (e) {
        if (!restored) {
          restored = true;
          await this.updateTree(snapshot).runInContext(ctx);
          this.events.log.next(LogEntry.error('Error during state transaction, reverting'));
        }
        if (isNested) {
          this.inTransactionError = true;
          throw e;
        }

        if (options?.rethrowErrors) {
          throw e;
        } else {
          console.error(e);
        }
      } finally {
        if (!isNested) {
          this.inTransaction = false;
          this.events.changed.next({ state: this, inTransaction: false });
          this.behaviors.isUpdating.next(false);

          if (!restored) {
            if (options?.canUndo)
              this.addHistory(snapshot, typeof options.canUndo === 'string' ? options.canUndo : void 0);
            else this.clearHistory();
          }
        }
      }
    });
  }

  private _inUpdate = false;
  /**
   * Determines whether the state is currently "inside" updateTree function.
   * This is different from "isUpdating" which wraps entire transactions.
   */
  get inUpdate() {
    return this._inUpdate;
  }

  /**
   * Queues up a reconciliation of the existing state tree.
   *
   * If the tree is StateBuilder.To<T>, the corresponding StateObject is returned by the task.
   * @param tree Tree instance or a tree builder instance
   * @param doNotReportTiming Indicates whether to log timing of the individual transforms
   */
  updateTree<T extends StateObject>(
    tree: StateBuilder.To<T, any>,
    options?: Partial<State.UpdateOptions>,
  ): Task<StateObjectSelector<T>>;
  updateTree(tree: StateTree | StateBuilder, options?: Partial<State.UpdateOptions>): Task<void>;
  updateTree(tree: StateTree | StateBuilder, options?: Partial<State.UpdateOptions>): Task<any> {
    const params: UpdateParams = { tree, options };
    return Task.create(
      'Update Tree',
      async (taskCtx) => {
        const removed = await this.updateQueue.enqueue(params);
        if (!removed) return;

        this._inUpdate = true;

        const snapshot = options?.canUndo ? this._tree.asImmutable() : void 0;
        let reverted = false;

        if (!this.inTransaction) this.behaviors.isUpdating.next(true);
        try {
          if (StateBuilder.is(tree)) {
            if (tree.editInfo.applied)
              throw new Error('This builder has already been applied. Create a new builder for further state updates');
            tree.editInfo.applied = true;
          }

          this.reverted = false;
          const ret =
            options && (options.revertIfAborted || options.revertOnError)
              ? await this._revertibleTreeUpdate(taskCtx, params, options)
              : await this._updateTree(taskCtx, params);
          reverted = this.reverted;

          if (ret.ctx.hadError) this.inTransactionError = true;

          if (!ret.cell) return;

          return new StateObjectSelector(ret.cell.transform.ref, this);
        } finally {
          this._inUpdate = false;
          this.updateQueue.handled(params);
          if (!this.inTransaction) {
            this.behaviors.isUpdating.next(false);
            if (!options?.canUndo) {
              if (!this.undoingHistory) this.clearHistory();
            } else if (!reverted) {
              this.addHistory(snapshot!, typeof options.canUndo === 'string' ? options.canUndo : void 0);
            }
          }
        }
      },
      () => {
        this.updateQueue.remove(params);
      },
    );
  }

  private reverted = false;
  private updateQueue = new AsyncQueue<UpdateParams>();
  private dependencyWaiters = new DependencyWaiters();

  private async _revertibleTreeUpdate(
    taskCtx: RuntimeContext,
    params: UpdateParams,
    options: Partial<State.UpdateOptions>,
  ) {
    const old = this.tree;
    const ret = await this._updateTree(taskCtx, params);
    const revert =
      ((ret.ctx.hadError || ret.ctx.wasAborted) && options.revertOnError) ||
      (ret.ctx.wasAborted && options.revertIfAborted);
    if (revert) {
      this.reverted = true;
      return await this._updateTree(taskCtx, { tree: old, options: params.options });
    }
    return ret;
  }

  private async _updateTree(taskCtx: RuntimeContext, params: UpdateParams) {
    let updated = false;
    const ctx = this.updateTreeAndCreateCtx(params.tree, taskCtx, params.options);
    try {
      updated = await update(ctx);
      if (StateBuilder.isTo(params.tree)) {
        const cell = this.select(params.tree.ref)[0];
        return { ctx, cell };
      }
      return { ctx };
    } finally {
      this.spine.current = undefined;

      if (updated) this.events.changed.next({ state: this, inTransaction: this.inTransaction });
    }
  }

  private updateTreeAndCreateCtx(
    tree: StateTree | StateBuilder,
    taskCtx: RuntimeContext,
    options: Partial<State.UpdateOptions> | undefined,
  ) {
    const _tree = (StateBuilder.is(tree) ? tree.getTree() : tree).asTransient();
    if (StateBuilder.is(tree)) this.carryFlagWrites(tree.editInfo, _tree);
    const oldTree = this._tree;
    this._tree = _tree;

    const ctx: UpdateContext = {
      parent: this,
      editInfo: StateBuilder.is(tree) ? tree.editInfo : void 0,

      errorFree: this.errorFree,
      dependencyWaiters: this.dependencyWaiters,
      flagLog: this.flagLog,
      orphans: void 0,
      taskCtx,
      oldTree,
      tree: _tree,
      cells: this.cells as Map<StateTransform.Ref, StateObjectCell>,
      spine: this.spine,

      results: [],

      options: { ...StateUpdateDefaultOptions, ...options },

      changed: false,
      hadError: false,
      wasAborted: false,
      newCurrent: void 0,

      getCellData: this.tryGetCellData,
    };

    this.errorFree = true;

    return ctx;
  }

  constructor(rootObject: StateObject, params: State.Params) {
    this._tree = StateTree.createEmpty(StateTransform.createRoot(params && params.rootState)).asTransient();
    const tree = this._tree;
    const root = tree.root;
    this.runTask = params.runTask;

    if (params?.historyCapacity !== void 0) this.historyCapacity = params.historyCapacity;

    (this.cells as Map<StateTransform.Ref, StateObjectCell>).set(root.ref, {
      parent: this,
      transform: root,
      sourceRef: void 0,
      obj: rootObject,
      status: 'ok',
      state: { ...root.state },
      errorText: void 0,
      params: {
        definition: {},
        values: {},
      },
      dependencies: { dependentBy: [], dependsOn: [] },
      cache: {},
    });

    this.globalContext = params && params.globalContext;
  }
}

namespace State {
  export interface Params {
    runTask<T>(task: Task<T>): Promise<T>;
    globalContext?: unknown;
    rootState?: StateTransform.State;
    historyCapacity?: number;
  }

  export function create(rootObject: StateObject, params: Params) {
    return new State(rootObject, params);
  }

  export type Cells = ReadonlyMap<StateTransform.Ref, StateObjectCell>;

  export type Tree = StateTree;
  export type Builder = StateBuilder;

  export interface ObjectEvent {
    state: State;
    ref: Ref;
  }

  export namespace ObjectEvent {
    export function isCell(e: ObjectEvent, cell?: StateObjectCell) {
      return !!cell && e.ref === cell.transform.ref && e.state === cell.parent;
    }
  }

  export interface Snapshot {
    readonly tree: StateTree.Serialized;
  }

  export interface UpdateOptions {
    doNotLogTiming: boolean;
    doNotUpdateCurrent: boolean;
    revertIfAborted: boolean;
    revertOnError: boolean;
    canUndo: boolean | string;
  }
}

const StateUpdateDefaultOptions: State.UpdateOptions = {
  doNotLogTiming: false,
  doNotUpdateCurrent: true,
  revertIfAborted: false,
  revertOnError: false,
  canUndo: false,
};

type Ref = StateTransform.Ref;

type UpdateParams = { tree: StateTree | StateBuilder; options?: Partial<State.UpdateOptions> };

interface UpdateContext {
  parent: State;
  editInfo: StateBuilder.EditInfo | undefined;

  errorFree: boolean;
  dependencyWaiters: DependencyWaiters;
  /** `State.flagLog`: deleted cells leave it. */
  flagLog: FlagLog;
  /** Cells whose dependency was deleted in this pass (`unlinkCell`): dependents, evaluated again. */
  orphans: StateObjectCell[] | undefined;
  taskCtx: RuntimeContext;
  oldTree: StateTree;
  tree: TransientTree;
  cells: Map<StateTransform.Ref, StateObjectCell>;
  spine: StateTreeSpine.Impl;

  results: UpdateNodeResult[];

  // suppress timing messages
  options: State.UpdateOptions;

  changed: boolean;
  hadError: boolean;
  wasAborted: boolean;
  newCurrent?: Ref;

  /**
   * Refs that are scheduled to be (re-)evaluated in this update pass:
   * the union of every `roots[i]` and its descendants, less those `doError`
   * cleared (below an error or a Null object). Used to distinguish
   * "dep not yet produced this pass" (defer + retry) from "dep can never
   * be produced this pass" (throw).
   */
  scheduled?: Set<Ref>;

  /**
   * Subtree roots that were skipped because at least one of their
   * dependencies hadn't been evaluated yet. Drained after the main loop
   * makes a fixpoint pass.
   */
  deferred?: Ref[];

  getCellData: (ref: string) => any;
}

async function update(ctx: UpdateContext) {
  // if only a single node was added/updated, we can skip potentially expensive diffing
  const fastTrack = !!(
    ctx.editInfo &&
    ctx.editInfo.count === 1 &&
    ctx.editInfo.lastUpdate &&
    ctx.editInfo.sourceTree === ctx.oldTree
  );
  // a builder on the current tree that moved no node, with the cells in line with that tree: the deletes of the
  // whole-tree diff are the subtrees it removed, the roots are among its edits (reconcile.md §2)
  const parent = ctx.parent as unknown as { cellsInSync: boolean };
  const cellsInSync = parent.cellsInSync;
  const edits = ctx.editInfo?.edits;
  const editsOnly = !fastTrack && cellsInSync && !!edits && !edits.moved && ctx.editInfo!.sourceTree === ctx.oldTree;
  // until this pass completes: a throw leaves cells behind the tree, which only the whole-tree diff finds
  parent.cellsInSync = false;

  let deletes: StateTransform.Ref[],
    deletedObjects: (StateObject | undefined)[] = [],
    roots: StateTransform.Ref[];

  if (fastTrack) {
    deletes = [];
    roots = [ctx.editInfo!.lastUpdate!];
  } else {
    // find all nodes that will definitely be deleted.
    // this is done in "post order", meaning that leaves will be deleted first.
    deletes = editsOnly
      ? findEditedDeletes(ctx.oldTree, ctx.tree, ctx.cells, edits!.removed)
      : findDeletes(ctx.oldTree, ctx.tree, ctx.cells);

    const current = ctx.parent.current;
    let hasCurrent = false;
    for (const d of deletes) {
      if (d === current) {
        hasCurrent = true;
        break;
      }
    }

    if (hasCurrent) {
      const newCurrent = findNewCurrent(ctx.oldTree, current, deletes, ctx.cells);
      ctx.parent.setCurrent(newCurrent);
    }

    for (let i = deletes.length - 1; i >= 0; i--) {
      const cell = ctx.cells.get(deletes[i]);
      if (cell) {
        dispose(cell.transform, cell.obj, cell?.transform.params, cell.cache, ctx.parent.globalContext);
      }
    }

    for (const d of deletes) {
      const cell = ctx.cells.get(d);
      if (cell) {
        cell.parent = void 0;
        unlinkCell(cell, ctx);
      }
      const obj = cell && cell.obj;
      ctx.cells.delete(d);
      ctx.flagLog.delete(d);
      deletedObjects.push(obj);
    }

    // Find roots where transform version changed or where nodes will be added.
    roots = editsOnly ? findEditedRoots(ctx.cells, ctx.tree, edits!.refs) : findUpdateRoots(ctx.cells, ctx.tree);
  }

  // Init empty cells where not present
  // this is done in "pre order", meaning that "parents" will be created 1st.
  const init = initCells(ctx, roots);

  // Notify additions of new cells.
  for (const cell of init.added) {
    ctx.parent.events.cell.created.next({ state: ctx.parent, ref: cell.transform.ref, cell });
  }

  for (let i = 0; i < deletes.length; i++) {
    const d = deletes[i];
    const parent = ctx.oldTree.transforms.get(d).parent;
    ctx.parent.events.object.removed.next({ state: ctx.parent, ref: d, obj: deletedObjects[i] });
    ctx.parent.events.cell.removed.next({ state: ctx.parent, ref: d, parent: parent });
  }

  if (deletedObjects.length) deletedObjects = [];

  if (init.dependent) {
    for (const cell of init.dependent) {
      roots.push(cell.transform.ref);
    }
  }

  // Set status of cells that will be updated to 'pending', and collect the refs that will be (re-)evaluated this
  // pass so that `updateSubtree` can distinguish "dep not produced yet" (defer + retry) from "dep can never be
  // produced" (throw via resolveDependencies). One walk for both (perf-audit-2 03/P4).
  initCellStatus(ctx, roots);

  // Sequentially update all the subtrees.
  for (const root of roots) {
    await updateSubtree(ctx, root);
  }

  // Drain the deferred queue: nodes whose `dependsOn` cells weren't yet
  // evaluated when first visited get retried until either all succeed or
  // a full pass makes no forward progress (true deadlock / unresolvable).
  if (ctx.deferred && ctx.deferred.length > 0) {
    while (ctx.deferred.length > 0) {
      const pending = ctx.deferred;
      ctx.deferred = [];
      let progress = false;
      for (const ref of pending) {
        const before = ctx.deferred.length;
        await updateSubtree(ctx, ref);
        // Forward progress = this attempt did not re-defer the same ref.
        const reDeferred = ctx.deferred.length > before && ctx.deferred[ctx.deferred.length - 1] === ref;
        if (!reDeferred) progress = true;
      }
      if (!progress) {
        const stuck = ctx.deferred.map((r) => {
          const c = ctx.cells.get(r);
          if (!c) return r;
          const blockers = pendingBlockers(ctx, c);
          return `${r} (waiting on: ${blockers.length ? blockers.join(', ') : 'unknown'})`;
        });
        ctx.deferred = [];
        throw new Error(`Unresolved dependency: ${stuck.join('; ')}`);
      }
    }
  }

  // Sync cell states
  if (!ctx.editInfo) {
    syncNewStates(ctx);
  }

  let newCurrent: StateTransform.Ref | undefined = ctx.newCurrent;
  // Raise object updated events
  for (const update of ctx.results) {
    if (update.action === 'created') {
      ctx.parent.events.object.created.next({ state: ctx.parent, ref: update.ref, obj: update.obj! });
      if (!ctx.newCurrent) {
        const transform = ctx.tree.transforms.get(update.ref);
        if (!transform.state.isGhost && update.obj !== StateObject.Null) newCurrent = update.ref;
      }
    } else if (update.action === 'updated') {
      ctx.parent.events.object.updated.next({
        state: ctx.parent,
        ref: update.ref,
        action: 'in-place',
        obj: update.obj,
        oldData: update.oldData,
      });
    } else if (update.action === 'replaced') {
      ctx.parent.events.object.updated.next({
        state: ctx.parent,
        ref: update.ref,
        action: 'recreate',
        obj: update.obj,
        oldObj: update.oldObj,
      });
    }
  }

  if (newCurrent) {
    if (!ctx.options.doNotUpdateCurrent) ctx.parent.setCurrent(newCurrent);
  } else {
    // check if old current or its parent hasn't become null
    const current = ctx.parent.current;
    const currentCell = ctx.cells.get(current);
    if (
      currentCell &&
      (currentCell.obj === StateObject.Null ||
        (currentCell.status === 'error' && currentCell.errorText === ParentNullErrorText))
    ) {
      newCurrent = findNewCurrent(ctx.oldTree, current, [], ctx.cells);
      ctx.parent.setCurrent(newCurrent);
    }
  }

  // the fast track evaluates one subtree: the cells stay in line only if the tree changed nowhere else (a write to the
  // builder's tree directly is not counted) and they were before
  parent.cellsInSync = fastTrack ? cellsInSync && editsWithin(ctx.tree, edits, ctx.editInfo!.lastUpdate!) : true;

  return deletes.length > 0 || roots.length > 0 || ctx.changed;
}

/** whether every recorded edit is in the subtree of `root` (none removed) */
function editsWithin(tree: StateTree, edits: TransientTree.Edits | undefined, root: Ref) {
  if (!edits || edits.removed.length > 0) return false;
  for (const ref of edits.refs) {
    let t = tree.transforms.get(ref);
    while (t && t.ref !== root && t.ref !== StateTransform.RootRef) t = tree.transforms.get(t.parent);
    if (!t || t.ref !== root) return false;
  }
  return true;
}

const StateFlags = ['isGhost', 'isLocked', 'isHidden', 'isCollapsed'] as const;
type FlagLog = Map<Ref, { [K in (typeof StateFlags)[number]]?: { n: number; value: boolean | undefined } }>;

function syncNewStatesVisitor(n: StateTransform, tree: StateTree, ctx: UpdateContext) {
  const cell = ctx.cells.get(n.ref);
  if (!cell || !StateTransform.syncState(cell.state, n.state)) return;
  ctx.parent.events.cell.stateUpdated.next({ state: ctx.parent, ref: n.ref, cell });
}

function syncNewStates(ctx: UpdateContext) {
  StateTree.doPreOrder(ctx.tree, ctx.tree.root, ctx, syncNewStatesVisitor);
}

function setCellStatus(ctx: UpdateContext, ref: Ref, status: StateObjectCell.Status, errorText?: string) {
  const cell = ctx.cells.get(ref)!;
  const changed = cell.status !== status;
  cell.status = status;
  cell.errorText = errorText;
  if (changed) ctx.parent.events.cell.stateUpdated.next({ state: ctx.parent, ref, cell });
}

function initCellStatusVisitor(t: StateTransform, _: any, ctx: UpdateContext) {
  const cell = ctx.cells.get(t.ref)!;
  cell.transform = t;
  if (relinkCells(cell, ctx)) {
    // Edges changed (e.g. param update added/removed dependencies) —
    // re-verify there is no cycle.
    checkDependenciesCycle(cell);
  }
  setCellStatus(ctx, t.ref, 'pending');
  ctx.scheduled!.add(t.ref);
}

function initCellStatus(ctx: UpdateContext, roots: Ref[]) {
  ctx.scheduled = new Set<Ref>();
  for (const root of roots) {
    StateTree.doPreOrder(ctx.tree, ctx.tree.transforms.get(root), ctx, initCellStatusVisitor);
  }
}

/** Removes a deleted cell from the dependency graph; the cells that depended on it wait for its ref again. */
function unlinkCell(cell: StateObjectCell, ctx: UpdateContext) {
  for (const other of cell.dependencies.dependsOn) {
    arraySetRemove(other.dependencies.dependentBy, cell);
  }
  const ref = cell.transform.ref;
  for (const by of cell.dependencies.dependentBy) {
    arraySetRemove(by.dependencies.dependsOn, cell);
    // a dependency derived from params is not removed with it (`TransientTree.remove` follows `dependsOn`)
    if (!by.parent) continue;
    ctx.dependencyWaiters.wait(by, ref);
    if (ctx.orphans) ctx.orphans.push(by);
    else ctx.orphans = [by];
  }
  ctx.dependencyWaiters.forget(cell);
}

type InitCellsCtx = { ctx: UpdateContext; visited: Set<Ref>; added: StateObjectCell[] };

function addCellsVisitor(transform: StateTransform, _: any, { ctx, added, visited }: InitCellsCtx) {
  visited.add(transform.ref);

  if (ctx.cells.has(transform.ref)) {
    return;
  }

  const cell: StateObjectCell = {
    parent: ctx.parent,
    transform,
    sourceRef: void 0,
    status: 'pending',
    state: { ...transform.state },
    errorText: void 0,
    params: void 0,
    dependencies: { dependentBy: [], dependsOn: [] },
    cache: void 0,
  };

  ctx.cells.set(transform.ref, cell);
  added.push(cell);
}

// type LinkCellsCtx = { ctx: UpdateContext, visited: Set<Ref>, dependent: UniqueArray<Ref, StateObjectCell> }

/**
 * Cells that depend on refs with no cell: a ref in the tree below a Null object or an error (nothing there is
 * evaluated), or a dependency deleted while its dependent stayed. `resolveDependencies` reports a waiting cell
 * unresolved, and `linkWaiting` links it once the ref has a cell (dependencies.md §3).
 */
class DependencyWaiters {
  private byRef = new Map<Ref, StateObjectCell[]>();
  private byCell = new Map<StateObjectCell, Ref[]>();

  isWaiting(cell: StateObjectCell) {
    return this.byCell.size > 0 && this.byCell.has(cell);
  }

  isWaitingFor(cell: StateObjectCell, ref: Ref) {
    return this.byCell.size > 0 && (this.byCell.get(cell)?.indexOf(ref) ?? -1) >= 0;
  }

  wait(cell: StateObjectCell, ref: Ref) {
    const cells = this.byRef.get(ref);
    if (cells) arraySetAdd(cells, cell);
    else this.byRef.set(ref, [cell]);
    const refs = this.byCell.get(cell);
    if (refs) arraySetAdd(refs, ref);
    else this.byCell.set(cell, [ref]);
  }

  /** The cells that waited for `ref`, which no longer wait for it. */
  take(ref: Ref): StateObjectCell[] | undefined {
    const cells = this.byRef.get(ref);
    if (!cells) return void 0;
    this.byRef.delete(ref);
    for (const cell of cells) this.removeRef(cell, ref);
    return cells;
  }

  /** Stop waiting for the refs `cell` no longer depends on. */
  retain(cell: StateObjectCell, refs: ReadonlySet<Ref>) {
    const waited = this.byCell.get(cell);
    if (!waited) return;
    for (let i = waited.length - 1; i >= 0; i--) {
      const ref = waited[i];
      if (refs.has(ref)) continue;
      waited.splice(i, 1);
      this.removeCell(ref, cell);
    }
    if (waited.length === 0) this.byCell.delete(cell);
  }

  forget(cell: StateObjectCell) {
    const waited = this.byCell.get(cell);
    if (!waited) return;
    this.byCell.delete(cell);
    for (const ref of waited) this.removeCell(ref, cell);
  }

  private removeRef(cell: StateObjectCell, ref: Ref) {
    const refs = this.byCell.get(cell);
    if (!refs) return;
    arraySetRemove(refs, ref);
    if (refs.length === 0) this.byCell.delete(cell);
  }

  private removeCell(ref: Ref, cell: StateObjectCell) {
    const cells = this.byRef.get(ref);
    if (!cells) return;
    arraySetRemove(cells, cell);
    if (cells.length === 0) this.byRef.delete(ref);
  }
}

/** Links `target` to the cell of `ref`, or makes it wait for one (`DependencyWaiters`). Returns `true` if it linked. */
function linkDependency(target: StateObjectCell, ref: Ref, ctx: UpdateContext) {
  const cell = ctx.cells.get(ref);
  if (!cell) {
    ctx.dependencyWaiters.wait(target, ref);
    return false;
  }
  arraySetAdd(target.dependencies.dependsOn, cell);
  arraySetAdd(cell.dependencies.dependentBy, target);
  return true;
}

/** Links the cells that waited for the new `cell` (`DependencyWaiters`); they are its dependents now. */
function linkWaiting(cell: StateObjectCell, ctx: UpdateContext, linked: StateObjectCell[]) {
  const waiting = ctx.dependencyWaiters.take(cell.transform.ref);
  if (!waiting) return;
  for (const w of waiting) {
    arraySetAdd(w.dependencies.dependsOn, cell);
    arraySetAdd(cell.dependencies.dependentBy, w);
    linked.push(w);
  }
}

function linkCells(target: StateObjectCell, ctx: UpdateContext) {
  const effective = StateTransform.getEffectiveDependsOn(target.transform);
  if (effective.length === 0) return;

  for (const ref of effective) {
    const t = ctx.tree.transforms.get(ref);
    if (!t) {
      throw new Error(`Cannot depend on a non-existent transform.`);
    }
    linkDependency(target, ref, ctx);
  }
}

/** Links newly added cells, and the cells that waited for them (returned), and checks for cycles. */
function linkAdded(ctx: UpdateContext, added: StateObjectCell[]) {
  const waited: StateObjectCell[] = [];
  for (const cell of added) {
    linkCells(cell, ctx);
  }
  for (const cell of added) {
    linkWaiting(cell, ctx, waited);
  }

  // Cycle detection over the dependency cross-edges of newly added cells.
  // Parent/child is already a tree, so cycles can only enter via dependsOn.
  for (const cell of added) {
    checkDependenciesCycle(cell);
  }
  for (const cell of waited) {
    checkDependenciesCycle(cell);
  }
  return waited;
}

/**
 * Diff the current outgoing dependency edges of `target` against the effective
 * set derived from its (possibly updated) transform. Unlinks stale edges and
 * links new ones so the dependency graph reflects param changes. Idempotent
 * when nothing has changed. Returns `true` if any edge was added or removed.
 */
function relinkCells(target: StateObjectCell, ctx: UpdateContext): boolean {
  const effective = StateTransform.getEffectiveDependsOn(target.transform);
  const current = target.dependencies.dependsOn;

  // Fast path: same number and all current refs are still in effective (and it waits for none).
  if (current.length === effective.length && !ctx.dependencyWaiters.isWaiting(target)) {
    // a set only pays off for many dependencies; most cells have a few
    const effectiveSet = effective.length > 8 ? new Set(effective) : undefined;
    let same = true;
    for (const c of current) {
      if (effectiveSet ? !effectiveSet.has(c.transform.ref) : effective.indexOf(c.transform.ref) < 0) {
        same = false;
        break;
      }
    }
    if (same) return false;
  }

  const desired = new Set(effective);
  let changed = false;
  ctx.dependencyWaiters.retain(target, desired);

  // Remove stale outgoing edges.
  for (let i = current.length - 1; i >= 0; i--) {
    const dep = current[i];
    if (!desired.has(dep.transform.ref)) {
      current.splice(i, 1);
      arraySetRemove(dep.dependencies.dependentBy, target);
      changed = true;
    }
  }

  // Add new outgoing edges.
  const have = new Set(current.map((c) => c.transform.ref));
  for (const ref of effective) {
    if (have.has(ref)) continue;
    const t = ctx.tree.transforms.get(ref);
    if (!t) {
      // a deleted dependency it waits for: still unresolved
      if (ctx.dependencyWaiters.isWaitingFor(target, ref)) continue;
      throw new Error(`Cannot depend on a non-existent transform.`);
    }
    if (linkDependency(target, ref, ctx)) changed = true;
  }

  return changed;
}

/**
 * Detect a cycle in the `dependsOn` graph reachable from `start`.
 *
 * Iterative DFS, returning the closed cycle path (first occurrence of the
 * repeated ref appended at the end) or `undefined` if none. Operates on the
 * already-linked cell graph; safe to call after `linkCells` / `relinkCells`.
 */
function detectDependenciesCycle(start: StateObjectCell): StateTransform.Ref[] | undefined {
  if (start.dependencies.dependsOn.length === 0) return void 0;
  const stack: { cell: StateObjectCell; idx: number }[] = [{ cell: start, idx: 0 }];
  const onPath = new Set<StateTransform.Ref>([start.transform.ref]);
  const fully = new Set<StateTransform.Ref>();

  while (stack.length > 0) {
    const top = stack[stack.length - 1];
    const deps = top.cell.dependencies.dependsOn;
    if (top.idx >= deps.length) {
      onPath.delete(top.cell.transform.ref);
      fully.add(top.cell.transform.ref);
      stack.pop();
      continue;
    }
    const next = deps[top.idx++];
    const ref = next.transform.ref;
    if (fully.has(ref)) continue;
    if (onPath.has(ref)) {
      const path: StateTransform.Ref[] = [];
      let started = false;
      for (const s of stack) {
        if (started || s.cell.transform.ref === ref) {
          started = true;
          path.push(s.cell.transform.ref);
        }
      }
      path.push(ref);
      return path;
    }
    onPath.add(ref);
    stack.push({ cell: next, idx: 0 });
  }
  return void 0;
}

function checkDependenciesCycle(cell: StateObjectCell) {
  const cycle = detectDependenciesCycle(cell);
  if (cycle) throw new StateTreeCycleError(cycle);
}

function initCells(ctx: UpdateContext, roots: Ref[]) {
  const initCtx: InitCellsCtx = { ctx, visited: new Set(), added: [] };

  // Add new cells
  for (const root of roots) {
    StateTree.doPreOrder(ctx.tree, ctx.tree.transforms.get(root), initCtx, addCellsVisitor);
  }

  linkAdded(ctx, initCtx.added);

  let dependent: UniqueArray<Ref, StateObjectCell>;

  // Find dependent cells
  initCtx.visited.forEach((ref) => {
    const cell = ctx.cells.get(ref)!;
    for (const by of cell.dependencies.dependentBy) {
      if (initCtx.visited.has(by.transform.ref)) continue;

      if (!dependent) dependent = UniqueArray.create();
      UniqueArray.add(dependent, by.transform.ref, by);
    }
  });

  // The cells whose dependency was deleted: their objects were made from it
  if (ctx.orphans) {
    for (const cell of ctx.orphans) {
      // deleted after it was left (later in the deletes), or in the subtrees above already
      if (ctx.cells.get(cell.transform.ref) !== cell || initCtx.visited.has(cell.transform.ref)) continue;
      if (!dependent!) dependent = UniqueArray.create();
      UniqueArray.add(dependent, cell.transform.ref, cell);
    }
  }

  // TODO: check if dependent cells are all "proper roots"

  // The subtrees of the dependents need cells too: a dependent below a Null object (where the diff does not descend)
  // can have children with none. A cell that waited for one of those is a dependent as well (reconcile.md §4).
  if (dependent!) {
    const depCtx: InitCellsCtx = { ctx, visited: initCtx.visited, added: [] };
    for (let i = 0; i < dependent.array.length; i++) {
      const start = depCtx.added.length;
      const root = ctx.tree.transforms.get(dependent.array[i].transform.ref);
      StateTree.doPreOrder(ctx.tree, root, depCtx, addCellsVisitor);
      for (const w of linkAdded(ctx, depCtx.added.slice(start))) {
        if (!depCtx.visited.has(w.transform.ref)) UniqueArray.add(dependent, w.transform.ref, w);
      }
    }
    for (const cell of depCtx.added) initCtx.added.push(cell);
  }

  return { added: initCtx.added, dependent: dependent! ? dependent!.array : void 0 };
}

function findNewCurrent(tree: StateTree, start: Ref, deletes: Ref[], cells: Map<Ref, StateObjectCell>) {
  const deleteSet = new Set(deletes);
  return _findNewCurrent(tree, start, deleteSet, cells);
}

function _findNewCurrent(tree: StateTree, ref: Ref, deletes: Set<Ref>, cells: Map<Ref, StateObjectCell>): Ref {
  if (ref === StateTransform.RootRef) return ref;

  const node = tree.transforms.get(ref)!;
  const siblings = tree.children.get(node.parent)!.values();

  let prevCandidate: Ref | undefined = void 0,
    seenRef = false;

  while (true) {
    const s = siblings.next();
    if (s.done) break;

    if (deletes.has(s.value)) continue;
    const cell = cells.get(s.value);
    if (!cell || cell.status === 'error' || cell.obj === StateObject.Null) {
      continue;
    }

    const t = tree.transforms.get(s.value);
    if (t.state.isGhost) continue;
    if (s.value === ref) {
      seenRef = true;
      if (!deletes.has(ref)) prevCandidate = ref;
      continue;
    }

    if (seenRef) return t.ref;

    prevCandidate = t.ref;
  }

  if (prevCandidate) return prevCandidate;
  return _findNewCurrent(tree, node.parent, deletes, cells);
}

/** Set status and error text of the cell. Remove all existing objects in the subtree. */
function doError(ctx: UpdateContext, ref: Ref, errorObject: any | undefined, silent: boolean) {
  if (!silent) {
    ctx.hadError = true;
    (ctx.parent as any as { errorFree: boolean }).errorFree = false;
  }

  const cell = ctx.cells.get(ref)!;
  // it will not produce an object this pass: what depends on it is not deferred for it (`pendingBlockers`)
  ctx.scheduled?.delete(ref);

  if (errorObject) {
    ctx.wasAborted = ctx.wasAborted || Task.isAbort(errorObject);
    const message = '' + errorObject;
    setCellStatus(ctx, ref, 'error', message);
    if (!silent) ctx.parent.events.log.next({ type: 'error', timestamp: new Date(), message });
  } else {
    cell.params = void 0;
  }

  if (cell.obj) {
    const obj = cell.obj;
    cell.obj = void 0;
    cell.cache = void 0;
    ctx.parent.events.object.removed.next({ state: ctx.parent, ref, obj });
  }

  // remove the objects in the child nodes if they exist
  const children = ctx.tree.children.get(ref).values();
  while (true) {
    const next = children.next();
    if (next.done) return;
    doError(ctx, next.value, void 0, silent);
  }
}

type UpdateNodeResult =
  | { ref: Ref; action: 'created'; obj: StateObject }
  | { ref: Ref; action: 'updated'; oldData?: any; obj: StateObject }
  | { ref: Ref; action: 'replaced'; oldObj?: StateObject; obj: StateObject }
  | { action: 'none' };

const ParentNullErrorText = 'Parent is null';

/**
 * Refs that `cell` depends on which are scheduled for evaluation in this
 * pass and not evaluated yet (still 'pending', with no object).
 * An empty result means deps are either ready, evaluated without an object
 * (in error, or below a Null object), or out-of-scope — in any of those cases
 * the cell should proceed (and may throw at `resolveDependencies` time if the
 * dep is genuinely missing). dependencies.md §4.
 */
function pendingBlockers(ctx: UpdateContext, cell: StateObjectCell): Ref[] {
  const blockers: Ref[] = [];
  const scheduled = ctx.scheduled;
  if (!scheduled) return blockers;
  for (const dep of cell.dependencies.dependsOn) {
    if (dep.obj) continue;
    // evaluated this pass and left without an object: waiting would not change that
    if (dep.status !== 'pending') continue;
    if (!scheduled.has(dep.transform.ref)) continue;
    blockers.push(dep.transform.ref);
  }
  return blockers;
}

async function updateSubtree(ctx: UpdateContext, root: Ref) {
  const cell = ctx.cells.get(root);
  if (cell && cell.dependencies.dependsOn.length > 0) {
    const blockers = pendingBlockers(ctx, cell);
    if (blockers.length > 0) {
      // A dependency hasn't been produced yet this pass — defer this
      // subtree and retry after other roots have run. Children will be
      // visited when the deferred re-attempt succeeds.
      if (!ctx.deferred) ctx.deferred = [];
      ctx.deferred.push(root);
      return;
    }
  }

  setCellStatus(ctx, root, 'processing');

  let isNull = false;
  try {
    const start = now();
    const update = await updateNode(ctx, root);
    const time = now() - start;

    if (update.action !== 'none') ctx.changed = true;

    setCellStatus(ctx, root, 'ok');
    ctx.results.push(update);
    if (update.action === 'created') {
      isNull = update.obj === StateObject.Null;
      if (!isNull && !ctx.options.doNotLogTiming)
        ctx.parent.events.log.next(LogEntry.info(`Created ${update.obj.label} in ${formatTimespan(time)}.`));
    } else if (update.action === 'updated') {
      isNull = update.obj === StateObject.Null;
      if (!isNull && !ctx.options.doNotLogTiming)
        ctx.parent.events.log.next(LogEntry.info(`Updated ${update.obj.label} in ${formatTimespan(time)}.`));
    } else if (update.action === 'replaced') {
      isNull = update.obj === StateObject.Null;
      if (!isNull && !ctx.options.doNotLogTiming)
        ctx.parent.events.log.next(LogEntry.info(`Updated ${update.obj.label} in ${formatTimespan(time)}.`));
    }
  } catch (e) {
    ctx.changed = true;
    if (!ctx.hadError) ctx.newCurrent = root;
    doError(ctx, root, e, false);
    console.error(e);
    return;
  }

  const children = ctx.tree.children.get(root).values();
  while (true) {
    const next = children.next();
    if (next.done) return;
    if (isNull) doError(ctx, next.value, void 0, true);
    else await updateSubtree(ctx, next.value);
  }
}

function resolveParams(ctx: UpdateContext, transform: StateTransform, src: StateObject, cell: StateObjectCell) {
  const prms = transform.transformer.definition.params;
  const definition = prms ? prms(src, ctx.parent.globalContext) : {};

  if (transform.version !== (transform as any)._normalized_param_version) {
    (transform.params as any) = ParamDefinition.normalizeParams(definition, transform.params, 'all');
    (transform as any)._normalized_param_version = transform.version;
  } else {
    const defaultValues = ParamDefinition.getDefaultValues(definition);
    (transform.params as any) = transform.params ? assignIfUndefined(transform.params, defaultValues) : defaultValues;
  }

  ParamDefinition.resolveRefs(definition, transform.params, ctx.getCellData);
  return { definition, values: transform.params };
}

async function updateNode(ctx: UpdateContext, currentRef: Ref): Promise<UpdateNodeResult> {
  const { oldTree, tree } = ctx;
  const current = ctx.cells.get(currentRef)!;
  const transform = current.transform;

  // Special case for Root
  if (current.transform.ref === StateTransform.RootRef) {
    return { action: 'none' };
  }

  const treeParent = ctx.cells.get(current.transform.parent);
  const isParentNull = treeParent?.obj === StateObject.Null;

  // Special case for when the immediate parent is null
  // This could happen then manually applying transforms to
  // already existing null nudes
  if (isParentNull) {
    current.sourceRef = treeParent.transform.ref;

    if (oldTree.transforms.has(currentRef) && current.params) {
      const oldParams = current.params.values;
      const oldCache = current.cache;
      dispose(transform, current.obj, oldParams, oldCache, ctx.parent.globalContext);

      current.params = undefined;
      current.obj = StateObject.Null;
      return { ref: currentRef, action: 'updated', obj: current.obj! };
    } else {
      current.params = undefined;
      return { ref: currentRef, action: 'created', obj: StateObject.Null };
    }
  }

  const parentCell =
    transform.transformer.definition.from.length === 0
      ? treeParent
      : StateSelection.findAncestorOfType(tree, ctx.cells, currentRef, transform.transformer.definition.from);
  if (!parentCell) {
    throw new Error(`No suitable parent found for '${currentRef}'`);
  }

  ctx.spine.current = current;

  const parent = parentCell.obj!;
  current.sourceRef = parentCell.transform.ref;

  const params = resolveParams(ctx, transform, parent, current);

  if (!oldTree.transforms.has(currentRef) || !current.params) {
    current.params = params;
    const obj = await createObject(ctx, current, transform.transformer, parent, params.values);
    updateTag(obj, transform);
    current.obj = obj;

    return { ref: currentRef, action: 'created', obj };
  } else {
    const oldParams = current.params.values;
    const oldCache = current.cache;
    const oldData = current.obj?.data;
    const newParams = params.values;
    current.params = params;

    const updateKind =
      !!current.obj && current.obj !== StateObject.Null
        ? await updateObject(ctx, current, transform.transformer, parent, current.obj!, oldParams, newParams)
        : StateTransformer.UpdateResult.Recreate;

    switch (updateKind) {
      case StateTransformer.UpdateResult.Recreate: {
        const oldObj = current.obj;
        dispose(transform, oldObj, oldParams, oldCache, ctx.parent.globalContext);

        const newObj = await createObject(ctx, current, transform.transformer, parent, newParams);

        updateTag(newObj, transform);
        current.obj = newObj;
        return { ref: currentRef, action: 'replaced', oldObj, obj: newObj };
      }
      case StateTransformer.UpdateResult.Updated:
        updateTag(current.obj, transform);
        return { ref: currentRef, action: 'updated', oldData, obj: current.obj! };
      case StateTransformer.UpdateResult.Null: {
        dispose(transform, current.obj, oldParams, oldCache, ctx.parent.globalContext);

        current.obj = StateObject.Null;
        return { ref: currentRef, action: 'updated', obj: current.obj! };
      }
      default:
        return { action: 'none' };
    }
  }
}

function dispose(transform: StateTransform, b: StateObject | undefined, params: any, cache: any, globalContext: any) {
  transform.transformer.definition.dispose?.(
    {
      b: b !== StateObject.Null ? b : void 0,
      params,
      cache,
    },
    globalContext,
  );
}

function updateTag(obj: StateObject | undefined, transform: StateTransform) {
  if (!obj || obj === StateObject.Null) return;
  (obj.tags as string[] | undefined) = transform.tags;
}

function runTask<T>(t: T | Task<T>, ctx: RuntimeContext) {
  if (typeof (t as any).runInContext === 'function') return (t as Task<T>).runInContext(ctx);
  return t as T;
}

function resolveDependencies(ctx: UpdateContext, cell: StateObjectCell) {
  // a dependency with no cell is not linked (`DependencyWaiters`)
  if (ctx.dependencyWaiters.isWaiting(cell)) throw new Error('Unresolved dependency.');
  if (cell.dependencies.dependsOn.length === 0) return void 0;

  const deps = Object.create(null);

  for (const dep of cell.dependencies.dependsOn) {
    if (!dep.obj) {
      throw new Error('Unresolved dependency.');
    }
    deps[dep.transform.ref] = dep.obj;
  }

  return deps;
}

function createObject(
  ctx: UpdateContext,
  cell: StateObjectCell,
  transformer: StateTransformer,
  a: StateObject,
  params: any,
) {
  if (!cell.cache) cell.cache = Object.create(null);
  return runTask(
    transformer.definition.apply(
      { a, params, cache: cell.cache, spine: ctx.spine, dependencies: resolveDependencies(ctx, cell) },
      ctx.parent.globalContext,
    ),
    ctx.taskCtx,
  );
}

async function updateObject(
  ctx: UpdateContext,
  cell: StateObjectCell,
  transformer: StateTransformer,
  a: StateObject,
  b: StateObject,
  oldParams: any,
  newParams: any,
) {
  if (!transformer.definition.update) {
    return StateTransformer.UpdateResult.Recreate;
  }
  if (!cell.cache) cell.cache = Object.create(null);
  return runTask(
    transformer.definition.update(
      { a, oldParams, b, newParams, cache: cell.cache, spine: ctx.spine, dependencies: resolveDependencies(ctx, cell) },
      ctx.parent.globalContext,
    ),
    ctx.taskCtx,
  );
}
