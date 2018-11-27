/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectCell } from './object';
import { StateTree } from './tree';
import { Transform } from './transform';
import { Transformer } from './transformer';
import { UUID } from 'mol-util';
import { RuntimeContext, Task } from 'mol-task';
import { StateSelection } from './state/selection';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { StateTreeBuilder } from './tree/builder';
import { StateAction } from './action';
import { StateActionManager } from './action/manager';
import { TransientTree } from './tree/transient';
import { LogEntry } from 'mol-util/log-entry';
import { now, formatTimespan } from 'mol-util/now';
import { ParamDefinition } from 'mol-util/param-definition';

export { State }

class State {
    private _tree: TransientTree = StateTree.createEmpty().asTransient();

    protected errorFree = true;
    private transformCache = new Map<Transform.Ref, unknown>();

    private ev = RxEventHelper.create();

    readonly globalContext: unknown = void 0;
    readonly events = {
        cell: {
            stateUpdated: this.ev<State.ObjectEvent & { cellState: StateObjectCell.State}>(),
            created: this.ev<State.ObjectEvent & { cell: StateObjectCell }>(),
            removed: this.ev<State.ObjectEvent & { parent: Transform.Ref }>(),
        },
        object: {
            updated: this.ev<State.ObjectEvent & { action: 'in-place' | 'recreate', obj: StateObject, oldObj?: StateObject }>(),
            created: this.ev<State.ObjectEvent & { obj: StateObject }>(),
            removed: this.ev<State.ObjectEvent & { obj?: StateObject }>()
        },
        log: this.ev<LogEntry>(),
        changed: this.ev<void>()
    };

    readonly behaviors = {
        currentObject: this.ev.behavior<State.ObjectEvent>({ state: this, ref: Transform.RootRef })
    };

    readonly actions = new StateActionManager();

    get tree(): StateTree { return this._tree; }
    get transforms() { return (this._tree as StateTree).transforms; }
    get cellStates() { return (this._tree as StateTree).cellStates; }
    get current() { return this.behaviors.currentObject.value.ref; }

    build() { return this._tree.build(); }

    readonly cells: State.Cells = new Map();

    getSnapshot(): State.Snapshot {
        return { tree: StateTree.toJSON(this._tree) };
    }

    setSnapshot(snapshot: State.Snapshot) {
        const tree = StateTree.fromJSON(snapshot.tree);
        return this.update(tree);
    }

    setCurrent(ref: Transform.Ref) {
        this.behaviors.currentObject.next({ state: this, ref });
    }

    updateCellState(ref: Transform.Ref, stateOrProvider: ((old: StateObjectCell.State) => Partial<StateObjectCell.State>) | Partial<StateObjectCell.State>) {
        const update = typeof stateOrProvider === 'function'
            ? stateOrProvider(this.tree.cellStates.get(ref))
            : stateOrProvider;

        if (this._tree.updateCellState(ref, update)) {
            this.events.cell.stateUpdated.next({ state: this, ref, cellState: this.tree.cellStates.get(ref) });
        }
    }

    dispose() {
        this.ev.dispose();
    }

    /**
     * Select Cells by ref or a query generated on the fly.
     * @example state.select('test')
     * @example state.select(q => q.byRef('test').subtree())
     */
    select(selector: Transform.Ref | ((q: typeof StateSelection.Generators) => StateSelection.Selector)) {
        if (typeof selector === 'string') return StateSelection.select(selector, this);
        return StateSelection.select(selector(StateSelection.Generators), this)
    }

    /** If no ref is specified, apply to root */
    apply<A extends StateAction>(action: A, params: StateAction.Params<A>, ref: Transform.Ref = Transform.RootRef): Task<void> {
        return Task.create('Apply Action', ctx => {
            const cell = this.cells.get(ref);
            if (!cell) throw new Error(`'${ref}' does not exist.`);
            if (cell.status !== 'ok') throw new Error(`Action cannot be applied to a cell with status '${cell.status}'`);

            return runTask(action.definition.run({ ref, cell, a: cell.obj!, params, state: this }, this.globalContext), ctx);
        });
    }

    update(tree: StateTree | StateTreeBuilder): Task<void> {
        const _tree = (StateTreeBuilder.is(tree) ? tree.getTree() : tree).asTransient();
        return Task.create('Update Tree', async taskCtx => {
            let updated = false;
            try {
                const oldTree = this._tree;
                this._tree = _tree;

                const ctx: UpdateContext = {
                    parent: this,
                    editInfo: StateTreeBuilder.is(tree) ? tree.editInfo : void 0,

                    errorFree: this.errorFree,
                    taskCtx,
                    oldTree,
                    tree: _tree,
                    cells: this.cells as Map<Transform.Ref, StateObjectCell>,
                    transformCache: this.transformCache,

                    results: [],

                    changed: false,
                    hadError: false,
                    newCurrent: void 0
                };

                this.errorFree = true;
                // TODO: handle "cancelled" error? Or would this be handled automatically?
                updated = await update(ctx);
            } finally {
                if (updated) this.events.changed.next();
            }
        });
    }

    constructor(rootObject: StateObject, params?: { globalContext?: unknown }) {
        const tree = this._tree;
        const root = tree.root;

        (this.cells as Map<Transform.Ref, StateObjectCell>).set(root.ref, {
            transform: root,
            sourceRef: void 0,
            obj: rootObject,
            status: 'ok',
            version: root.version,
            errorText: void 0
        });

        this.globalContext = params && params.globalContext;
    }
}

namespace State {
    export type Cells = ReadonlyMap<Transform.Ref, StateObjectCell>

    export type Tree = StateTree
    export type Builder = StateTreeBuilder

    export interface ObjectEvent {
        state: State,
        ref: Ref
    }

    export interface Snapshot {
        readonly tree: StateTree.Serialized
    }

    export function create(rootObject: StateObject, params?: { globalContext?: unknown, defaultObjectProps?: unknown }) {
        return new State(rootObject, params);
    }
}

type Ref = Transform.Ref

interface UpdateContext {
    parent: State,
    editInfo: StateTreeBuilder.EditInfo | undefined

    errorFree: boolean,
    taskCtx: RuntimeContext,
    oldTree: StateTree,
    tree: TransientTree,
    cells: Map<Transform.Ref, StateObjectCell>,
    transformCache: Map<Ref, unknown>,

    results: UpdateNodeResult[],

    changed: boolean,
    hadError: boolean,
    newCurrent?: Ref
}

async function update(ctx: UpdateContext) {
    // if only a single node was added/updated, we can skip potentially expensive diffing
    const fastTrack = !!(ctx.errorFree && ctx.editInfo && ctx.editInfo.count === 1 && ctx.editInfo.lastUpdate && ctx.editInfo.sourceTree === ctx.oldTree);

    let deletes: Transform.Ref[], deletedObjects: (StateObject | undefined)[] = [], roots: Transform.Ref[];

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

        for (const d of deletes) {
            const obj = ctx.cells.has(d) ? ctx.cells.get(d)!.obj : void 0;
            ctx.cells.delete(d);
            ctx.transformCache.delete(d);
            deletedObjects.push(obj);
        }

        // Find roots where transform version changed or where nodes will be added.
        roots = findUpdateRoots(ctx.cells, ctx.tree);
    }

    // Init empty cells where not present
    // this is done in "pre order", meaning that "parents" will be created 1st.
    const addedCells = initCells(ctx, roots);

    // Ensure cell states stay consistent
    if (!ctx.editInfo) {
        syncStates(ctx);
    }

    // Notify additions of new cells.
    for (const cell of addedCells) {
        ctx.parent.events.cell.created.next({ state: ctx.parent, ref: cell.transform.ref, cell });
    }

    for (let i = 0; i < deletes.length; i++) {
        const d = deletes[i];
        const parent = ctx.oldTree.transforms.get(d).parent;
        ctx.parent.events.object.removed.next({ state: ctx.parent, ref: d, obj: deletedObjects[i] });
        ctx.parent.events.cell.removed.next({ state: ctx.parent, ref: d, parent: parent });
    }

    if (deletedObjects.length) deletedObjects = [];

    // Set status of cells that will be updated to 'pending'.
    initCellStatus(ctx, roots);

    // Sequentially update all the subtrees.
    for (const root of roots) {
        await updateSubtree(ctx, root);
    }

    let newCurrent: Transform.Ref | undefined = ctx.newCurrent;
    // Raise object updated events
    for (const update of ctx.results) {
        if (update.action === 'created') {
            ctx.parent.events.object.created.next({ state: ctx.parent, ref: update.ref, obj: update.obj! });
            if (!ctx.newCurrent) {
                const transform = ctx.tree.transforms.get(update.ref);
                if (!(transform.props && transform.props.isGhost) && update.obj !== StateObject.Null) newCurrent = update.ref;
            }
        } else if (update.action === 'updated') {
            ctx.parent.events.object.updated.next({ state: ctx.parent, ref: update.ref, action: 'in-place', obj: update.obj });
        } else if (update.action === 'replaced') {
            ctx.parent.events.object.updated.next({ state: ctx.parent, ref: update.ref, action: 'recreate', obj: update.obj, oldObj: update.oldObj });
        }
    }

    if (newCurrent) ctx.parent.setCurrent(newCurrent);
    else {
        // check if old current or its parent hasn't become null
        const current = ctx.parent.current;
        const currentCell = ctx.cells.get(current);
        if (currentCell && (
                currentCell.obj === StateObject.Null
            || (currentCell.status === 'error' && currentCell.errorText === ParentNullErrorText))) {
            newCurrent = findNewCurrent(ctx.oldTree, current, [], ctx.cells);
            ctx.parent.setCurrent(newCurrent);
        }
    }


    return deletes.length > 0 || roots.length > 0 || ctx.changed;
}

function findUpdateRoots(cells: Map<Transform.Ref, StateObjectCell>, tree: StateTree) {
    const findState = { roots: [] as Ref[], cells };
    StateTree.doPreOrder(tree, tree.root, findState, findUpdateRootsVisitor);
    return findState.roots;
}

function findUpdateRootsVisitor(n: Transform, _: any, s: { roots: Ref[], cells: Map<Ref, StateObjectCell> }) {
    const cell = s.cells.get(n.ref);
    if (!cell || cell.version !== n.version || cell.status === 'error') {
        s.roots.push(n.ref);
        return false;
    }
    // nothing below a Null object can be an update root
    if (cell && cell.obj === StateObject.Null) return false;
    return true;
}

type FindDeletesCtx = { newTree: StateTree, cells: State.Cells, deletes: Ref[] }
function checkDeleteVisitor(n: Transform, _: any, ctx: FindDeletesCtx) {
    if (!ctx.newTree.transforms.has(n.ref) && ctx.cells.has(n.ref)) ctx.deletes.push(n.ref);
}
function findDeletes(ctx: UpdateContext): Ref[] {
    const deleteCtx: FindDeletesCtx = { newTree: ctx.tree, cells: ctx.cells, deletes: [] };
    StateTree.doPostOrder(ctx.oldTree, ctx.oldTree.root, deleteCtx, checkDeleteVisitor);
    return deleteCtx.deletes;
}

function syncStatesVisitor(n: Transform, tree: StateTree, oldState: StateTree.CellStates) {
    if (!oldState.has(n.ref)) return;
    (tree as TransientTree).updateCellState(n.ref, oldState.get(n.ref));
}
function syncStates(ctx: UpdateContext) {
    StateTree.doPreOrder(ctx.tree, ctx.tree.root, ctx.oldTree.cellStates, syncStatesVisitor);
}

function setCellStatus(ctx: UpdateContext, ref: Ref, status: StateObjectCell.Status, errorText?: string) {
    const cell = ctx.cells.get(ref)!;
    const changed = cell.status !== status;
    cell.status = status;
    cell.errorText = errorText;
    if (changed) ctx.parent.events.cell.stateUpdated.next({ state: ctx.parent, ref, cellState: ctx.tree.cellStates.get(ref) });
}

function initCellStatusVisitor(t: Transform, _: any, ctx: UpdateContext) {
    ctx.cells.get(t.ref)!.transform = t;
    setCellStatus(ctx, t.ref, 'pending');
}

function initCellStatus(ctx: UpdateContext, roots: Ref[]) {
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.transforms.get(root), ctx, initCellStatusVisitor);
    }
}

type InitCellsCtx = { ctx: UpdateContext, added: StateObjectCell[] }
function initCellsVisitor(transform: Transform, _: any, { ctx, added }: InitCellsCtx) {
    if (ctx.cells.has(transform.ref)) {
        return;
    }

    const cell: StateObjectCell = {
        transform,
        sourceRef: void 0,
        status: 'pending',
        version: UUID.create22(),
        errorText: void 0
    };
    ctx.cells.set(transform.ref, cell);
    added.push(cell);
}

function initCells(ctx: UpdateContext, roots: Ref[]) {
    const initCtx: InitCellsCtx = { ctx, added: [] };
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.transforms.get(root), initCtx, initCellsVisitor);
    }
    return initCtx.added;
}

function findNewCurrent(tree: StateTree, start: Ref, deletes: Ref[], cells: Map<Ref, StateObjectCell>) {
    const deleteSet = new Set(deletes);
    return _findNewCurrent(tree, start, deleteSet, cells);
}

function _findNewCurrent(tree: StateTree, ref: Ref, deletes: Set<Ref>, cells: Map<Ref, StateObjectCell>): Ref {
    if (ref === Transform.RootRef) return ref;

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
        if (t.props && t.props.isGhost) continue;
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
function doError(ctx: UpdateContext, ref: Ref, errorText: string | undefined, silent: boolean) {
    if (!silent) {
        ctx.hadError = true;
        (ctx.parent as any as { errorFree: boolean }).errorFree = false;
    }

    if (errorText) {
        setCellStatus(ctx, ref, 'error', errorText);
        if (!silent) ctx.parent.events.log.next({ type: 'error', timestamp: new Date(), message: errorText });
    }

    const cell = ctx.cells.get(ref)!;
    if (cell.obj) {
        const obj = cell.obj;
        cell.obj = void 0;
        ctx.parent.events.object.removed.next({ state: ctx.parent, ref, obj });
        ctx.transformCache.delete(ref);
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
            if (!isNull) ctx.parent.events.log.next(LogEntry.info(`Created ${update.obj.label} in ${formatTimespan(time)}.`));
        } else if (update.action === 'updated') {
            isNull = update.obj === StateObject.Null;
            if (!isNull) ctx.parent.events.log.next(LogEntry.info(`Updated ${update.obj.label} in ${formatTimespan(time)}.`));
        } else if (update.action === 'replaced') {
            isNull = update.obj === StateObject.Null;
            if (!isNull) ctx.parent.events.log.next(LogEntry.info(`Updated ${update.obj.label} in ${formatTimespan(time)}.`));
        }
    } catch (e) {
        ctx.changed = true;
        if (!ctx.hadError) ctx.newCurrent = root;
        doError(ctx, root, '' + e, false);
        return;
    }

    const children = ctx.tree.children.get(root).values();
    while (true) {
        const next = children.next();
        if (next.done) return;
        if (isNull) doError(ctx, next.value, ParentNullErrorText, true);
        else await updateSubtree(ctx, next.value);
    }
}

function resolveDefaultParams(ctx: UpdateContext, transform: Transform, src: StateObject) {
    const prms = transform.transformer.definition.params;
    const defaults = prms
        ? ParamDefinition.getDefaultValues(prms(src, ctx.parent.globalContext))
        : { };
    // TODO: this should probably be resolved each time separately the transform is applied.
    // the params should be cached in the cell?
    (transform.params as any) = defaults;
}

async function updateNode(ctx: UpdateContext, currentRef: Ref): Promise<UpdateNodeResult> {
    const { oldTree, tree } = ctx;
    const current = ctx.cells.get(currentRef)!;
    const transform = current.transform;

    // special case for Root
    if (current.transform.ref === Transform.RootRef) return { action: 'none' };

    const parentCell = StateSelection.findAncestorOfType(tree, ctx.cells, currentRef, transform.transformer.definition.from);
    if (!parentCell) {
        throw new Error(`No suitable parent found for '${currentRef}'`);
    }

    const parent = parentCell.obj!;
    current.sourceRef = parentCell.transform.ref;

    if (!transform.params) {
        resolveDefaultParams(ctx, transform, parent);
    }

    if (!oldTree.transforms.has(currentRef)) {
        const obj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
        current.obj = obj;
        current.version = transform.version;

        return { ref: currentRef, action: 'created', obj };
    } else {
        const oldParams = oldTree.transforms.get(currentRef)!.params;

        const updateKind = !!current.obj && current.obj !== StateObject.Null
            ? await updateObject(ctx, currentRef, transform.transformer, parent, current.obj!, oldParams, transform.params)
            : Transformer.UpdateResult.Recreate;

        switch (updateKind) {
            case Transformer.UpdateResult.Recreate: {
                const oldObj = current.obj;
                const newObj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
                current.obj = newObj;
                current.version = transform.version;
                return { ref: currentRef, action: 'replaced', oldObj, obj: newObj };
            }
            case Transformer.UpdateResult.Updated:
                current.version = transform.version;
                return { ref: currentRef, action: 'updated', obj: current.obj! };
            default:
                return { action: 'none' };
        }
    }
}

function runTask<T>(t: T | Task<T>, ctx: RuntimeContext) {
    if (typeof (t as any).runInContext === 'function') return (t as Task<T>).runInContext(ctx);
    return t as T;
}

function createObject(ctx: UpdateContext, ref: Ref, transformer: Transformer, a: StateObject, params: any) {
    const cache = Object.create(null);
    ctx.transformCache.set(ref, cache);
    return runTask(transformer.definition.apply({ a, params, cache }, ctx.parent.globalContext), ctx.taskCtx);
}

async function updateObject(ctx: UpdateContext, ref: Ref, transformer: Transformer, a: StateObject, b: StateObject, oldParams: any, newParams: any) {
    if (!transformer.definition.update) {
        return Transformer.UpdateResult.Recreate;
    }
    let cache = ctx.transformCache.get(ref);
    if (!cache) {
        cache = Object.create(null);
        ctx.transformCache.set(ref, cache);
    }
    return runTask(transformer.definition.update({ a, oldParams, b, newParams, cache }, ctx.parent.globalContext), ctx.taskCtx);
}