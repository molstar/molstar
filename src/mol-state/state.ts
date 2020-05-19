/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectCell, StateObjectSelector } from './object';
import { StateTree } from './tree';
import { StateTransform } from './transform';
import { StateTransformer } from './transformer';
import { RuntimeContext, Task } from '../mol-task';
import { StateSelection } from './state/selection';
import { RxEventHelper } from '../mol-util/rx-event-helper';
import { StateBuilder } from './state/builder';
import { StateAction } from './action';
import { StateActionManager } from './action/manager';
import { TransientTree } from './tree/transient';
import { LogEntry } from '../mol-util/log-entry';
import { now, formatTimespan } from '../mol-util/now';
import { ParamDefinition } from '../mol-util/param-definition';
import { StateTreeSpine } from './tree/spine';
import { AsyncQueue } from '../mol-util/async-queue';
import { isProductionMode } from '../mol-util/debug';
import { arraySetAdd, arraySetRemove } from '../mol-util/array';
import { UniqueArray } from '../mol-data/generic';
import { assignIfUndefined } from '../mol-util/object';

export { State };

class State {
    private _tree: TransientTree;

    protected errorFree = true;

    private ev = RxEventHelper.create();

    readonly globalContext: unknown = void 0;
    readonly events = {
        cell: {
            stateUpdated: this.ev<State.ObjectEvent & { cell: StateObjectCell }>(),
            created: this.ev<State.ObjectEvent & { cell: StateObjectCell }>(),
            removed: this.ev<State.ObjectEvent & { parent: StateTransform.Ref }>(),
        },
        object: {
            updated: this.ev<State.ObjectEvent & { action: 'in-place' | 'recreate', obj: StateObject, oldObj?: StateObject }>(),
            created: this.ev<State.ObjectEvent & { obj: StateObject }>(),
            removed: this.ev<State.ObjectEvent & { obj?: StateObject }>()
        },
        log: this.ev<LogEntry>(),
        changed: this.ev<{ state: State, inTransaction: boolean }>(),
        historyUpdated: this.ev<{ state: State }>()
    };

    readonly behaviors = {
        currentObject: this.ev.behavior<State.ObjectEvent>({ state: this, ref: StateTransform.RootRef }),
        isUpdating: this.ev.behavior<boolean>(false),
    };

    readonly actions = new StateActionManager();

    readonly runTask: <T>(task: Task<T>) => Promise<T>;

    get tree(): StateTree { return this._tree; }
    get transforms() { return (this._tree as StateTree).transforms; }
    get current() { return this.behaviors.currentObject.value.ref; }
    get root() { return this.cells.get((this._tree as StateTree).root.ref)!; }

    build() { return new StateBuilder.Root(this.tree, this); }

    readonly cells: State.Cells = new Map();
    private spine = new StateTreeSpine.Impl(this.cells);

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
        return Task.create('Undo', async ctx => {
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

    setSnapshot(snapshot: State.Snapshot) {
        const tree = StateTree.fromJSON(snapshot.tree);
        return this.updateTree(tree);
    }

    setCurrent(ref: StateTransform.Ref) {
        this.behaviors.currentObject.next({ state: this, ref });
    }

    updateCellState(ref: StateTransform.Ref, stateOrProvider: ((old: StateTransform.State) => Partial<StateTransform.State>) | Partial<StateTransform.State>) {
        const cell = this.cells.get(ref);
        if (!cell) return;

        const update = typeof stateOrProvider === 'function' ? stateOrProvider(cell.state) : stateOrProvider;

        if (StateTransform.assignState(cell.state, update)) {
            cell.transform = this._tree.assignState(cell.transform.ref, update);
            this.events.cell.stateUpdated.next({ state: this, ref, cell });
        }
    }

    dispose() {
        this.ev.dispose();
        this.actions.dispose();
    }

    /**
     * Select Cells using the provided selector.
     * @example state.query(StateSelection.Generators.byRef('test').ancestorOfType([type]))
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
    applyAction<A extends StateAction>(action: A, params: StateAction.Params<A>, ref: StateTransform.Ref = StateTransform.RootRef): Task<void> {
        return Task.create('Apply Action', ctx => {
            const cell = this.cells.get(ref);
            if (!cell) throw new Error(`'${ref}' does not exist.`);
            if (cell.status !== 'ok') throw new Error(`Action cannot be applied to a cell with status '${cell.status}'`);

            return runTask(action.definition.run({ ref, cell, a: cell.obj!, params, state: this }, this.globalContext), ctx);
        });
    }

    private inTransaction = false;
    private inTransactionError = false;

    /** Apply series of updates to the state. If any of them fail, revert to the original state. */
    transaction(edits: (ctx: RuntimeContext) => Promise<void> | void, options?: { canUndo?: string | boolean }) {
        return Task.create('State Transaction', async ctx => {
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
                    this.events.log.error(e);
                }
                if (isNested) {
                    this.inTransactionError = true;
                    throw e;
                }
            } finally {
                if (!isNested) {
                    this.inTransaction = false;
                    this.events.changed.next({ state: this, inTransaction: false });
                    this.behaviors.isUpdating.next(false);

                    if (!restored) {
                        if (options?.canUndo) this.addHistory(snapshot, typeof options.canUndo === 'string' ? options.canUndo : void 0);
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
    get inUpdate() { return this._inUpdate; }

    /**
     * Queues up a reconciliation of the existing state tree.
     *
     * If the tree is StateBuilder.To<T>, the corresponding StateObject is returned by the task.
     * @param tree Tree instance or a tree builder instance
     * @param doNotReportTiming Indicates whether to log timing of the individual transforms
     */
    updateTree<T extends StateObject>(tree: StateBuilder.To<T, any>, options?: Partial<State.UpdateOptions>): Task<StateObjectSelector<T>>
    updateTree(tree: StateTree | StateBuilder, options?: Partial<State.UpdateOptions>): Task<void>
    updateTree(tree: StateTree | StateBuilder, options?: Partial<State.UpdateOptions>): Task<any> {
        const params: UpdateParams = { tree, options };
        return Task.create('Update Tree', async taskCtx => {
            const removed = await this.updateQueue.enqueue(params);
            if (!removed) return;

            this._inUpdate = true;

            const snapshot = options?.canUndo ? this._tree.asImmutable() : void 0;
            let reverted = false;

            if (!this.inTransaction) this.behaviors.isUpdating.next(true);
            try {
                if (StateBuilder.is(tree)) {
                    if (tree.editInfo.applied) throw new Error('This builder has already been applied. Create a new builder for further state updates');
                    tree.editInfo.applied = true;
                }

                this.reverted = false;
                const ret = options && (options.revertIfAborted || options.revertOnError)
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
        }, () => {
            this.updateQueue.remove(params);
        });
    }

    private reverted = false;
    private updateQueue = new AsyncQueue<UpdateParams>();

    private async _revertibleTreeUpdate(taskCtx: RuntimeContext, params: UpdateParams, options: Partial<State.UpdateOptions>) {
        const old = this.tree;
        const ret = await this._updateTree(taskCtx, params);
        let revert = ((ret.ctx.hadError || ret.ctx.wasAborted) && options.revertOnError) || (ret.ctx.wasAborted && options.revertIfAborted);
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

    private updateTreeAndCreateCtx(tree: StateTree | StateBuilder, taskCtx: RuntimeContext, options: Partial<State.UpdateOptions> | undefined) {
        const _tree = (StateBuilder.is(tree) ? tree.getTree() : tree).asTransient();
        const oldTree = this._tree;
        this._tree = _tree;

        const ctx: UpdateContext = {
            parent: this,
            editInfo: StateBuilder.is(tree) ? tree.editInfo : void 0,

            errorFree: this.errorFree,
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
            newCurrent: void 0
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
                values: {}
            },
            dependencies: { dependentBy: [], dependsOn: [] },
            cache: { }
        });

        this.globalContext = params && params.globalContext;
    }
}

namespace State {
    export interface Params {
        runTask<T>(task: Task<T>): Promise<T>,
        globalContext?: unknown,
        rootState?: StateTransform.State,
        historyCapacity?: number
    }

    export function create(rootObject: StateObject, params: Params) {
        return new State(rootObject, params);
    }

    export type Cells = ReadonlyMap<StateTransform.Ref, StateObjectCell>

    export type Tree = StateTree
    export type Builder = StateBuilder

    export interface ObjectEvent {
        state: State,
        ref: Ref
    }

    export namespace ObjectEvent {
        export function isCell(e: ObjectEvent, cell?: StateObjectCell) {
            return !!cell && e.ref === cell.transform.ref && e.state === cell.parent;
        }
    }

    export interface Snapshot {
        readonly tree: StateTree.Serialized
    }

    export interface UpdateOptions {
        doNotLogTiming: boolean,
        doNotUpdateCurrent: boolean,
        revertIfAborted: boolean,
        revertOnError: boolean,
        canUndo: boolean | string
    }
}

const StateUpdateDefaultOptions: State.UpdateOptions = {
    doNotLogTiming: false,
    doNotUpdateCurrent: true,
    revertIfAborted: false,
    revertOnError: false,
    canUndo: false
};

type Ref = StateTransform.Ref

type UpdateParams = { tree: StateTree | StateBuilder, options?: Partial<State.UpdateOptions> }

interface UpdateContext {
    parent: State,
    editInfo: StateBuilder.EditInfo | undefined

    errorFree: boolean,
    taskCtx: RuntimeContext,
    oldTree: StateTree,
    tree: TransientTree,
    cells: Map<StateTransform.Ref, StateObjectCell>,
    spine: StateTreeSpine.Impl,

    results: UpdateNodeResult[],

    // suppress timing messages
    options: State.UpdateOptions,

    changed: boolean,
    hadError: boolean,
    wasAborted: boolean,
    newCurrent?: Ref
}

async function update(ctx: UpdateContext) {
    // if only a single node was added/updated, we can skip potentially expensive diffing
    const fastTrack = !!(ctx.editInfo && ctx.editInfo.count === 1 && ctx.editInfo.lastUpdate && ctx.editInfo.sourceTree === ctx.oldTree);
    let deletes: StateTransform.Ref[], deletedObjects: (StateObject | undefined)[] = [], roots: StateTransform.Ref[];

    if (fastTrack) {
        deletes = [];
        roots = [ctx.editInfo!.lastUpdate!];
    } else {
        // find all nodes that will definitely be deleted.
        // this is done in "post order", meaning that leaves will be deleted first.
        deletes = findDeletes(ctx);

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
                unlinkCell(cell);
            }
            const obj = cell && cell.obj;
            ctx.cells.delete(d);
            deletedObjects.push(obj);
        }

        // Find roots where transform version changed or where nodes will be added.
        roots = findUpdateRoots(ctx.cells, ctx.tree);
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

    // Set status of cells that will be updated to 'pending'.
    initCellStatus(ctx, roots);

    // Sequentially update all the subtrees.
    for (const root of roots) {
        await updateSubtree(ctx, root);
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
            ctx.parent.events.object.updated.next({ state: ctx.parent, ref: update.ref, action: 'in-place', obj: update.obj });
        } else if (update.action === 'replaced') {
            ctx.parent.events.object.updated.next({ state: ctx.parent, ref: update.ref, action: 'recreate', obj: update.obj, oldObj: update.oldObj });
        }
    }

    if (newCurrent) {
        if (!ctx.options.doNotUpdateCurrent) ctx.parent.setCurrent(newCurrent);
    } else {
        // check if old current or its parent hasn't become null
        const current = ctx.parent.current;
        const currentCell = ctx.cells.get(current);
        if (currentCell && (currentCell.obj === StateObject.Null
            || (currentCell.status === 'error' && currentCell.errorText === ParentNullErrorText))) {
            newCurrent = findNewCurrent(ctx.oldTree, current, [], ctx.cells);
            ctx.parent.setCurrent(newCurrent);
        }
    }

    return deletes.length > 0 || roots.length > 0 || ctx.changed;
}

function findUpdateRoots(cells: Map<StateTransform.Ref, StateObjectCell>, tree: StateTree) {
    const findState = { roots: [] as Ref[], cells };
    StateTree.doPreOrder(tree, tree.root, findState, findUpdateRootsVisitor);
    return findState.roots;
}

function findUpdateRootsVisitor(n: StateTransform, _: any, s: { roots: Ref[], cells: Map<Ref, StateObjectCell> }) {
    const cell = s.cells.get(n.ref);
    if (!cell || cell.transform.version !== n.version) {
        s.roots.push(n.ref);
        return false;
    }
    if (cell.status === 'error') return false;

    // nothing below a Null object can be an update root
    if (cell && cell.obj === StateObject.Null) return false;
    return true;
}

type FindDeletesCtx = { newTree: StateTree, cells: State.Cells, deletes: Ref[] }
function checkDeleteVisitor(n: StateTransform, _: any, ctx: FindDeletesCtx) {
    if (!ctx.newTree.transforms.has(n.ref) && ctx.cells.has(n.ref)) ctx.deletes.push(n.ref);
}
function findDeletes(ctx: UpdateContext): Ref[] {
    const deleteCtx: FindDeletesCtx = { newTree: ctx.tree, cells: ctx.cells, deletes: [] };
    StateTree.doPostOrder(ctx.oldTree, ctx.oldTree.root, deleteCtx, checkDeleteVisitor);
    return deleteCtx.deletes;
}

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
    ctx.cells.get(t.ref)!.transform = t;
    setCellStatus(ctx, t.ref, 'pending');
}

function initCellStatus(ctx: UpdateContext, roots: Ref[]) {
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.transforms.get(root), ctx, initCellStatusVisitor);
    }
}

function unlinkCell(cell: StateObjectCell) {
    for (const other of cell.dependencies.dependsOn) {
        arraySetRemove(other.dependencies.dependentBy, cell);
    }
}

type InitCellsCtx = { ctx: UpdateContext, visited: Set<Ref>, added: StateObjectCell[] }

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
        cache: void 0
    };

    ctx.cells.set(transform.ref, cell);
    added.push(cell);
}

// type LinkCellsCtx = { ctx: UpdateContext, visited: Set<Ref>, dependent: UniqueArray<Ref, StateObjectCell> }

function linkCells(target: StateObjectCell, ctx: UpdateContext) {
    if (!target.transform.dependsOn) return;

    for (const ref of target.transform.dependsOn) {
        const t = ctx.tree.transforms.get(ref);
        if (!t) {
            throw new Error(`Cannot depend on a non-existent transform.`);
        }

        const cell = ctx.cells.get(ref)!;
        arraySetAdd(target.dependencies.dependsOn, cell);
        arraySetAdd(cell.dependencies.dependentBy, target);
    }
}

function initCells(ctx: UpdateContext, roots: Ref[]) {
    const initCtx: InitCellsCtx = { ctx, visited: new Set(), added: [] };

    // Add new cells
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.transforms.get(root), initCtx, addCellsVisitor);
    }

    // Update links for newly added cells
    for (const cell of initCtx.added) {
        linkCells(cell, ctx);
    }

    let dependent: UniqueArray<Ref, StateObjectCell>;

    // Find dependent cells
    initCtx.visited.forEach(ref => {
        const cell = ctx.cells.get(ref)!;
        for (const by of cell.dependencies.dependentBy) {
            if (initCtx.visited.has(by.transform.ref)) continue;

            if (!dependent) dependent = UniqueArray.create();
            UniqueArray.add(dependent, by.transform.ref, by);
        }
    });

    // TODO: check if dependent cells are all "proper roots"

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

    let prevCandidate: Ref | undefined = void 0, seenRef = false;

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
    | { ref: Ref, action: 'created', obj: StateObject }
    | { ref: Ref, action: 'updated', obj: StateObject }
    | { ref: Ref, action: 'replaced', oldObj?: StateObject, obj: StateObject }
    | { action: 'none' }

const ParentNullErrorText = 'Parent is null';

async function updateSubtree(ctx: UpdateContext, root: Ref) {
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
            if (!isNull && !ctx.options.doNotLogTiming) ctx.parent.events.log.next(LogEntry.info(`Created ${update.obj.label} in ${formatTimespan(time)}.`));
        } else if (update.action === 'updated') {
            isNull = update.obj === StateObject.Null;
            if (!isNull && !ctx.options.doNotLogTiming) ctx.parent.events.log.next(LogEntry.info(`Updated ${update.obj.label} in ${formatTimespan(time)}.`));
        } else if (update.action === 'replaced') {
            isNull = update.obj === StateObject.Null;
            if (!isNull && !ctx.options.doNotLogTiming) ctx.parent.events.log.next(LogEntry.info(`Updated ${update.obj.label} in ${formatTimespan(time)}.`));
        }
    } catch (e) {
        ctx.changed = true;
        if (!ctx.hadError) ctx.newCurrent = root;
        doError(ctx, root, e, false);
        if (!isProductionMode) console.error(e);
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

function resolveParams(ctx: UpdateContext, transform: StateTransform, src: StateObject) {
    const prms = transform.transformer.definition.params;
    const definition = prms ? prms(src, ctx.parent.globalContext) : {};
    const defaultValues = ParamDefinition.getDefaultValues(definition);
    (transform.params as any) = transform.params
        ? assignIfUndefined(transform.params, defaultValues)
        : defaultValues;
    return { definition, values: transform.params };
}

async function updateNode(ctx: UpdateContext, currentRef: Ref): Promise<UpdateNodeResult> {
    const { oldTree, tree } = ctx;
    const current = ctx.cells.get(currentRef)!;
    const transform = current.transform;

    // special case for Root
    if (current.transform.ref === StateTransform.RootRef) {
        return { action: 'none' };
    }

    let parentCell = transform.transformer.definition.from.length === 0
        ? ctx.cells.get(current.transform.parent)
        : StateSelection.findAncestorOfType(tree, ctx.cells, currentRef, transform.transformer.definition.from);
    if (!parentCell) {
        throw new Error(`No suitable parent found for '${currentRef}'`);
    }

    ctx.spine.current = current;

    const parent = parentCell.obj!;
    current.sourceRef = parentCell.transform.ref;

    const params = resolveParams(ctx, transform, parent);

    if (!oldTree.transforms.has(currentRef) || !current.params) {
        current.params = params;
        const obj = await createObject(ctx, current, transform.transformer, parent, params.values);
        updateTag(obj, transform);
        current.obj = obj;

        return { ref: currentRef, action: 'created', obj };
    } else {
        const oldParams = current.params.values;
        const oldCache = current.cache;
        const newParams = params.values;
        current.params = params;


        const updateKind = !!current.obj && current.obj !== StateObject.Null
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
                return { ref: currentRef, action: 'updated', obj: current.obj! };
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
    transform.transformer.definition.dispose?.({
        b: b !== StateObject.Null ? b : void 0,
        params,
        cache
    }, globalContext);
}

function updateTag(obj: StateObject | undefined, transform: StateTransform) {
    if (!obj || obj === StateObject.Null) return;
    (obj.tags as string[] | undefined) = transform.tags;
}

function runTask<T>(t: T | Task<T>, ctx: RuntimeContext) {
    if (typeof (t as any).runInContext === 'function') return (t as Task<T>).runInContext(ctx);
    return t as T;
}

function resolveDependencies(cell: StateObjectCell) {
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

function createObject(ctx: UpdateContext, cell: StateObjectCell, transformer: StateTransformer, a: StateObject, params: any) {
    if (!cell.cache) cell.cache = Object.create(null);
    return runTask(transformer.definition.apply({ a, params, cache: cell.cache, spine: ctx.spine, dependencies: resolveDependencies(cell) }, ctx.parent.globalContext), ctx.taskCtx);
}

async function updateObject(ctx: UpdateContext, cell: StateObjectCell,  transformer: StateTransformer, a: StateObject, b: StateObject, oldParams: any, newParams: any) {
    if (!transformer.definition.update) {
        return StateTransformer.UpdateResult.Recreate;
    }
    if (!cell.cache) cell.cache = Object.create(null);
    return runTask(transformer.definition.update({ a, oldParams, b, newParams, cache: cell.cache, spine: ctx.spine, dependencies: resolveDependencies(cell) }, ctx.parent.globalContext), ctx.taskCtx);
}