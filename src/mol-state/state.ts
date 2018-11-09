/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectCell } from './object';
import { StateTree } from './tree';
import { Transform } from './transform';
import { Transformer } from './transformer';
import { StateContext } from './context';
import { UUID } from 'mol-util';
import { RuntimeContext, Task } from 'mol-task';
import { StateSelection } from './selection';

export { State }

class State {
    private _tree: StateTree = StateTree.createEmpty();
    private _current: Transform.Ref = this._tree.root.ref;
    private transformCache = new Map<Transform.Ref, unknown>();

    get tree() { return this._tree; }
    get current() { return this._current; }

    readonly cells: State.Cells = new Map();
    readonly context: StateContext;

    getSnapshot(): State.Snapshot {
        return { tree: StateTree.toJSON(this._tree) };
    }

    setSnapshot(snapshot: State.Snapshot) {
        const tree = StateTree.fromJSON(snapshot.tree);
        return this.update(tree);
    }

    setCurrent(ref: Transform.Ref) {
        this._current = ref;
        this.context.behaviors.currentObject.next({ ref });
    }

    updateCellState(ref: Transform.Ref, state?: Partial<StateObjectCell.State>) {
        // TODO
    }

    dispose() {
        this.context.dispose();
    }

    select(selector: StateSelection.Selector) {
        return StateSelection.select(selector, this);
    }

    get selector() {
        return StateSelection;
    }

    update(tree: StateTree): Task<void> {
        // TODO: support props
        return Task.create('Update Tree', async taskCtx => {
            try {
                const oldTree = this._tree;
                this._tree = tree;

                const ctx: UpdateContext = {
                    stateCtx: this.context,
                    taskCtx,
                    oldTree,
                    tree,
                    cells: this.cells,
                    transformCache: this.transformCache
                };
                // TODO: have "cancelled" error? Or would this be handled automatically?
                await update(ctx);
            } finally {
                this.context.events.updated.next();
            }
        });
    }

    constructor(rootObject: StateObject, params?: { globalContext?: unknown }) {
        const tree = this._tree;
        const root = tree.root;

        this.cells.set(root.ref, {
            ref: root.ref,
            obj: rootObject,
            status: 'ok',
            version: root.version,
            state: { ...StateObjectCell.DefaultState }
        });

        this.context = new StateContext({
            globalContext: params && params.globalContext
        });
    }
}

namespace State {
    export type Cells = Map<Transform.Ref, StateObjectCell>

    export interface Snapshot {
        readonly tree: StateTree.Serialized
    }

    export function create(rootObject: StateObject, params?: { globalContext?: unknown, defaultObjectProps?: unknown }) {
        return new State(rootObject, params);
    }
}

type Ref = Transform.Ref

interface UpdateContext {
    stateCtx: StateContext,
    taskCtx: RuntimeContext,
    oldTree: StateTree,
    tree: StateTree,
    cells: State.Cells,
    transformCache: Map<Ref, unknown>
}

async function update(ctx: UpdateContext) {
    const roots = findUpdateRoots(ctx.cells, ctx.tree);
    const deletes = findDeletes(ctx);
    for (const d of deletes) {
        const obj = ctx.cells.has(d) ? ctx.cells.get(d)!.obj : void 0;
        ctx.cells.delete(d);
        ctx.transformCache.delete(d);
        ctx.stateCtx.events.object.removed.next({ ref: d, obj });
        // TODO: handle current object change
    }

    initCells(ctx, roots);
    initObjectStatus(ctx, roots);

    for (const root of roots) {
        await updateSubtree(ctx, root);
    }
}

function findUpdateRoots(cells: State.Cells, tree: StateTree) {
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

function setObjectStatus(ctx: UpdateContext, ref: Ref, status: StateObjectCell.Status, errorText?: string) {
    const obj = ctx.cells.get(ref)!;
    const changed = obj.status !== status;
    obj.status = status;
    obj.errorText = errorText;
    if (changed) ctx.stateCtx.events.object.stateChanged.next({ ref });
}

function _initObjectStatusVisitor(t: Transform, _: any, ctx: UpdateContext) {
    setObjectStatus(ctx, t.ref, 'pending');
}

/** Return "resolve set" */
function initObjectStatus(ctx: UpdateContext, roots: Ref[]) {
    for (const root of roots) {
        StateTree.doPreOrder(ctx.tree, ctx.tree.nodes.get(root), ctx, _initObjectStatusVisitor);
    }
}

function _initCellsVisitor(transform: Transform, _: any, ctx: UpdateContext) {
    if (ctx.cells.has(transform.ref)) return;

    const obj: StateObjectCell = {
        ref: transform.ref,
        status: 'pending',
        version: UUID.create(),
        errorText: void 0,
        state: { ...StateObjectCell.DefaultState, ...transform.cellState }
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
    setObjectStatus(ctx, ref, 'error', errorText);
    const wrap = ctx.cells.get(ref)!;
    if (wrap.obj) {
        ctx.stateCtx.events.object.removed.next({ ref });
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

function findAncestor(tree: StateTree, cells: State.Cells, root: Ref, types: { type: StateObject.Type }[]): StateObject {
    let current = tree.nodes.get(root)!;
    while (true) {
        current = tree.nodes.get(current.parent)!;
        if (current.ref === Transform.RootRef) {
            return cells.get(Transform.RootRef)!.obj!;
        }
        const obj = cells.get(current.ref)!.obj!;
        for (const t of types) if (obj.type === t.type) return cells.get(current.ref)!.obj!;
    }
}

async function updateSubtree(ctx: UpdateContext, root: Ref) {
    setObjectStatus(ctx, root, 'processing');

    try {
        const update = await updateNode(ctx, root);
        setObjectStatus(ctx, root, 'ok');
        if (update.action === 'created') {
            ctx.stateCtx.events.object.created.next({ ref: root, obj: update.obj! });
        } else if (update.action === 'updated') {
            ctx.stateCtx.events.object.updated.next({ ref: root, obj: update.obj });
        } else if (update.action === 'replaced') {
            ctx.stateCtx.events.object.replaced.next({ ref: root, oldObj: update.oldObj, newObj: update.newObj });
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
    const { oldTree, tree, cells } = ctx;
    const transform = tree.nodes.get(currentRef);
    const parent = findAncestor(tree, cells, currentRef, transform.transformer.definition.from);
    // console.log('parent', transform.transformer.id, transform.transformer.definition.from[0].type, parent ? parent.ref : 'undefined')
    if (!oldTree.nodes.has(currentRef) || !cells.has(currentRef)) {
        // console.log('creating...', transform.transformer.id, oldTree.nodes.has(currentRef), objects.has(currentRef));
        const obj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
        const cell = cells.get(currentRef)!;
        cell.obj = obj;
        cell.version = transform.version;

        return { action: 'created', obj };
    } else {
        const current = cells.get(currentRef)!;
        const oldParams = oldTree.nodes.get(currentRef)!.params;

        const updateKind = current.status === 'ok' || current.ref === Transform.RootRef
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
    return runTask(transformer.definition.apply({ a, params, cache }, ctx.stateCtx.globalContext), ctx.taskCtx);
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
    return runTask(transformer.definition.update({ a, oldParams, b, newParams, cache }, ctx.stateCtx.globalContext), ctx.taskCtx);
}