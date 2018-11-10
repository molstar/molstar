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

export { State }

class State {
    private _tree: StateTree = StateTree.createEmpty();
    private _current: Transform.Ref = this._tree.root.ref;
    private transformCache = new Map<Transform.Ref, unknown>();

    private ev = RxEventHelper.create();

    readonly globalContext: unknown = void 0;
    readonly events = {
        object: {
            cellState: this.ev<State.ObjectEvent & { cell: StateObjectCell }>(),

            updated: this.ev<State.ObjectEvent & { obj?: StateObject }>(),
            replaced: this.ev<State.ObjectEvent & { oldObj?: StateObject, newObj?: StateObject }>(),
            created: this.ev<State.ObjectEvent & { obj: StateObject }>(),
            removed: this.ev<State.ObjectEvent & { obj?: StateObject }>()
        },
        warn: this.ev<string>(),
        updated: this.ev<void>()
    };

    readonly behaviors = {
        currentObject: this.ev.behavior<State.ObjectEvent>({ state: this, ref: Transform.RootRef })
    };

    get tree() { return this._tree; }
    get current() { return this._current; }

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
        this._current = ref;
        this.behaviors.currentObject.next({ state: this, ref });
    }

    updateCellState(ref: Transform.Ref, state?: Partial<StateObjectCell.State>) {
        // TODO
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

    apply(action: StateAction, ref: Transform.Ref) {
        
    }

    update(tree: StateTree | StateTreeBuilder): Task<void> {
        // TODO: support cell state
        const _tree = StateTreeBuilder.is(tree) ? tree.getTree() : tree;
        return Task.create('Update Tree', async taskCtx => {
            try {
                const oldTree = this._tree;
                this._tree = _tree;

                const ctx: UpdateContext = {
                    parent: this,
                    taskCtx,
                    oldTree,
                    tree: _tree,
                    cells: this.cells as Map<Transform.Ref, StateObjectCell>,
                    transformCache: this.transformCache
                };
                // TODO: have "cancelled" error? Or would this be handled automatically?
                await update(ctx);
            } finally {
                this.events.updated.next();
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
            version: root.version
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

    taskCtx: RuntimeContext,
    oldTree: StateTree,
    tree: StateTree,
    cells: Map<Transform.Ref, StateObjectCell>,
    transformCache: Map<Ref, unknown>
}

async function update(ctx: UpdateContext) {
    const roots = findUpdateRoots(ctx.cells, ctx.tree);
    const deletes = findDeletes(ctx);
    for (const d of deletes) {
        const obj = ctx.cells.has(d) ? ctx.cells.get(d)!.obj : void 0;
        ctx.cells.delete(d);
        ctx.transformCache.delete(d);
        ctx.parent.events.object.removed.next({ state: ctx.parent, ref: d, obj });
        // TODO: handle current object change
    }

    initCells(ctx, roots);
    initCellStatus(ctx, roots);

    for (const root of roots) {
        await updateSubtree(ctx, root);
    }
}

function findUpdateRoots(cells: Map<Transform.Ref, StateObjectCell>, tree: StateTree) {
    const findState = { roots: [] as Ref[], cells };
    StateTree.doPreOrder(tree, tree.root, findState, _findUpdateRoots);
    return findState.roots;
}

function _findUpdateRoots(n: Transform, _: any, s: { roots: Ref[], cells: Map<Ref, StateObjectCell> }) {
    const cell = s.cells.get(n.ref);
    if (!cell || cell.version !== n.version) {
        s.roots.push(n.ref);
        return false;
    }
    return true;
}

type FindDeletesCtx = { newTree: StateTree, cells: State.Cells, deletes: Ref[] }
function _visitCheckDelete(n: Transform, _: any, ctx: FindDeletesCtx) {
    if (!ctx.newTree.nodes.has(n.ref) && ctx.cells.has(n.ref)) ctx.deletes.push(n.ref);
}
function findDeletes(ctx: UpdateContext): Ref[] {
    const deleteCtx: FindDeletesCtx = { newTree: ctx.tree, cells: ctx.cells, deletes: [] };
    StateTree.doPostOrder(ctx.oldTree, ctx.oldTree.root, deleteCtx, _visitCheckDelete);
    return deleteCtx.deletes;
}

function setCellStatus(ctx: UpdateContext, ref: Ref, status: StateObjectCell.Status, errorText?: string) {
    const cell = ctx.cells.get(ref)!;
    const changed = cell.status !== status;
    cell.status = status;
    cell.errorText = errorText;
    if (changed) ctx.parent.events.object.cellState.next({ state: ctx.parent, ref, cell });
}

function _initCellStatusVisitor(t: Transform, _: any, ctx: UpdateContext) {
    ctx.cells.get(t.ref)!.transform = t;
    setCellStatus(ctx, t.ref, 'pending');
}

/** Return "resolve set" */
function initCellStatus(ctx: UpdateContext, roots: Ref[]) {
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.nodes.get(root), ctx, _initCellStatusVisitor);
    }
}

function _initCellsVisitor(transform: Transform, _: any, ctx: UpdateContext) {
    if (ctx.cells.has(transform.ref)) return;

    const obj: StateObjectCell = {
        transform,
        sourceRef: void 0,
        status: 'pending',
        version: UUID.create(),
        errorText: void 0
    };
    ctx.cells.set(transform.ref, obj);

    // TODO: created event???
}

function initCells(ctx: UpdateContext, roots: Ref[]) {
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.nodes.get(root), ctx, _initCellsVisitor);
    }
}

function doError(ctx: UpdateContext, ref: Ref, errorText: string) {
    setCellStatus(ctx, ref, 'error', errorText);
    const wrap = ctx.cells.get(ref)!;
    if (wrap.obj) {
        ctx.parent.events.object.removed.next({ state: ctx.parent, ref });
        ctx.transformCache.delete(ref);
        wrap.obj = void 0;
    }

    const children = ctx.tree.children.get(ref).values();
    while (true) {
        const next = children.next();
        if (next.done) return;
        doError(ctx, next.value, 'Parent node contains error.');
    }
}

async function updateSubtree(ctx: UpdateContext, root: Ref) {
    setCellStatus(ctx, root, 'processing');

    try {
        const update = await updateNode(ctx, root);
        setCellStatus(ctx, root, 'ok');
        if (update.action === 'created') {
            ctx.parent.events.object.created.next({ state: ctx.parent, ref: root, obj: update.obj! });
        } else if (update.action === 'updated') {
            ctx.parent.events.object.updated.next({ state: ctx.parent, ref: root, obj: update.obj });
        } else if (update.action === 'replaced') {
            ctx.parent.events.object.replaced.next({ state: ctx.parent, ref: root, oldObj: update.oldObj, newObj: update.newObj });
        }
    } catch (e) {
        doError(ctx, root, '' + e);
        return;
    }

    const children = ctx.tree.children.get(root).values();
    while (true) {
        const next = children.next();
        if (next.done) return;
        await updateSubtree(ctx, next.value);
    }
}

async function updateNode(ctx: UpdateContext, currentRef: Ref) {
    const { oldTree, tree } = ctx;
    const transform = tree.nodes.get(currentRef);
    const parentCell = StateSelection.findAncestorOfType(tree, ctx.cells, currentRef, transform.transformer.definition.from);

    if (!parentCell) {
        throw new Error(`No suitable parent found for '${currentRef}'`);
    }

    const parent = parentCell.obj!;
    const current = ctx.cells.get(currentRef)!;
    current.sourceRef = parentCell.transform.ref;

    // console.log('parent', transform.transformer.id, transform.transformer.definition.from[0].type, parent ? parent.ref : 'undefined')
    if (!oldTree.nodes.has(currentRef)) {
        // console.log('creating...', transform.transformer.id, oldTree.nodes.has(currentRef), objects.has(currentRef));
        const obj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
        current.obj = obj;
        current.version = transform.version;

        return { action: 'created', obj };
    } else {
        const oldParams = oldTree.nodes.get(currentRef)!.params;

        const updateKind = current.status === 'ok' || current.transform.ref === Transform.RootRef
            ? await updateObject(ctx, currentRef, transform.transformer, parent, current.obj!, oldParams, transform.params)
            : Transformer.UpdateResult.Recreate;

        switch (updateKind) {
            case Transformer.UpdateResult.Recreate: {
                const oldObj = current.obj;
                const newObj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
                current.obj = newObj;
                current.version = transform.version;
                return { action: 'replaced', oldObj, newObj: newObj };
            }
            case Transformer.UpdateResult.Updated:
                current.version = transform.version;
                return { action: 'updated', obj: current.obj };
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
    const cache = {};
    ctx.transformCache.set(ref, cache);
    return runTask(transformer.definition.apply({ a, params, cache }, ctx.parent.globalContext), ctx.taskCtx);
}

async function updateObject(ctx: UpdateContext, ref: Ref, transformer: Transformer, a: StateObject, b: StateObject, oldParams: any, newParams: any) {
    if (!transformer.definition.update) {
        return Transformer.UpdateResult.Recreate;
    }
    let cache = ctx.transformCache.get(ref);
    if (!cache) {
        cache = {};
        ctx.transformCache.set(ref, cache);
    }
    return runTask(transformer.definition.update({ a, oldParams, b, newParams, cache }, ctx.parent.globalContext), ctx.taskCtx);
}