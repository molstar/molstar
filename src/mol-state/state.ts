/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject, StateObjectCell } from './object';
import { StateTree } from './tree';
import { Transform } from './transform';
import { ImmutableTree } from './immutable-tree';
import { Transformer } from './transformer';
import { StateContext } from './context';
import { UUID } from 'mol-util';
import { RuntimeContext, Task } from 'mol-task';

export { State }

class State {
    private _tree: StateTree = StateTree.create();
    private _current: Transform.Ref = this._tree.rootRef;
    private transformCache = new Map<Transform.Ref, unknown>();

    get tree() { return this._tree; }
    get current() { return this._current; }

    readonly cells: State.Cells = new Map();
    readonly context: StateContext;

    getSnapshot(): State.Snapshot {
        const props = Object.create(null);
        const keys = this.cells.keys();
        while (true) {
            const key = keys.next();
            if (key.done) break;
            const o = this.cells.get(key.value)!;
            props[key.value] = { ...o.state };
        }
        return {
            tree: StateTree.toJSON(this._tree),
            props
        };
    }

    setSnapshot(snapshot: State.Snapshot) {
        const tree = StateTree.fromJSON(snapshot.tree);
        // TODO: support props and async
        return this.update(tree).run();
    }

    setCurrent(ref: Transform.Ref) {
        this._current = ref;
        this.context.behaviors.currentObject.next({ ref });
    }

    dispose() {
        this.context.dispose();
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
                    tree: tree,
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

    constructor(rootObject: StateObject, params?: { globalContext?: unknown, defaultCellState?: unknown }) {
        const tree = this._tree;
        const root = tree.getValue(tree.rootRef)!;
        const defaultCellState = (params && params.defaultCellState) || { }

        this.cells.set(tree.rootRef, {
            ref: tree.rootRef,
            obj: rootObject,
            status: 'ok',
            version: root.version,
            state: { ...defaultCellState }
        });

        this.context = new StateContext({
            globalContext: params && params.globalContext,
            defaultCellState,
            rootRef: tree.rootRef
        });
    }
}

namespace State {
    export type Cells = Map<Transform.Ref, StateObjectCell>

    export interface Snapshot {
        readonly tree: StateTree.Serialized,
        readonly props: { [key: string]: unknown }
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

        initObjectState(ctx, roots);

        for (const root of roots) {
            await updateSubtree(ctx, root);
        }
    }

    function findUpdateRoots(objects: State.Cells, tree: StateTree) {
        const findState = {
            roots: [] as Ref[],
            objects
        };

        ImmutableTree.doPreOrder(tree, tree.nodes.get(tree.rootRef)!, findState, (n, _, s) => {
            if (!s.objects.has(n.ref)) {
                s.roots.push(n.ref);
                return false;
            }
            const o = s.objects.get(n.ref)!;
            if (o.version !== n.value.version) {
                s.roots.push(n.ref);
                return false;
            }

            return true;
        });

        return findState.roots;
    }

    type FindDeletesCtx = { newTree: StateTree, cells: State.Cells, deletes: Ref[] }
    function _visitCheckDelete(n: ImmutableTree.Node<any>, _: any, ctx: FindDeletesCtx) {
        if (!ctx.newTree.nodes.has(n.ref) && ctx.cells.has(n.ref)) ctx.deletes.push(n.ref);
    }
    function findDeletes(ctx: UpdateContext): Ref[] {
        const deleteCtx: FindDeletesCtx = { newTree: ctx.tree, cells: ctx.cells, deletes: [] };
        ImmutableTree.doPostOrder(ctx.oldTree, ctx.oldTree.nodes.get(ctx.oldTree.rootRef), deleteCtx, _visitCheckDelete);
        return deleteCtx.deletes;
    }

    function setObjectState(ctx: UpdateContext, ref: Ref, status: StateObjectCell.Status, errorText?: string) {
        let changed = false;
        if (ctx.cells.has(ref)) {
            const obj = ctx.cells.get(ref)!;
            changed = obj.status !== status;
            obj.status = status;
            obj.errorText = errorText;
        } else {
            const obj: StateObjectCell = { ref, status, version: UUID.create(), errorText, state: { ...ctx.stateCtx.defaultCellState } };
            ctx.cells.set(ref, obj);
            changed = true;
        }
        if (changed) ctx.stateCtx.events.object.stateChanged.next({ ref });
    }

    function _initVisitor(t: ImmutableTree.Node<Transform>, _: any, ctx: UpdateContext) {
        setObjectState(ctx, t.ref, 'pending');
    }
    /** Return "resolve set" */
    function initObjectState(ctx: UpdateContext, roots: Ref[]) {
        for (const root of roots) {
            ImmutableTree.doPreOrder(ctx.tree, ctx.tree.nodes.get(root), ctx, _initVisitor);
        }
    }

    function doError(ctx: UpdateContext, ref: Ref, errorText: string) {
        setObjectState(ctx, ref, 'error', errorText);
        const wrap = ctx.cells.get(ref)!;
        if (wrap.obj) {
            ctx.stateCtx.events.object.removed.next({ ref });
            ctx.transformCache.delete(ref);
            wrap.obj = void 0;
        }

        const children = ctx.tree.nodes.get(ref)!.children.values();
        while (true) {
            const next = children.next();
            if (next.done) return;
            doError(ctx, next.value, 'Parent node contains error.');
        }
    }

    function findAncestor(tree: StateTree, objects: State.Cells, root: Ref, types: { type: StateObject.Type }[]): StateObject {
        let current = tree.nodes.get(root)!;
        while (true) {
            current = tree.nodes.get(current.parent)!;
            if (current.ref === tree.rootRef) {
                return objects.get(tree.rootRef)!.obj!;
            }
            const obj = objects.get(current.ref)!.obj!;
            for (const t of types) if (obj.type === t.type) return objects.get(current.ref)!.obj!;
        }
    }

    async function updateSubtree(ctx: UpdateContext, root: Ref) {
        setObjectState(ctx, root, 'processing');

        try {
            const update = await updateNode(ctx, root);
            setObjectState(ctx, root, 'ok');
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

        const children = ctx.tree.nodes.get(root)!.children.values();
        while (true) {
            const next = children.next();
            if (next.done) return;
            await updateSubtree(ctx, next.value);
        }
    }

    async function updateNode(ctx: UpdateContext, currentRef: Ref) {
        const { oldTree, tree, cells } = ctx;
        const transform = tree.getValue(currentRef)!;
        const parent = findAncestor(tree, cells, currentRef, transform.transformer.definition.from);
        // console.log('parent', transform.transformer.id, transform.transformer.definition.from[0].type, parent ? parent.ref : 'undefined')
        if (!oldTree.nodes.has(currentRef) || !cells.has(currentRef)) {
            // console.log('creating...', transform.transformer.id, oldTree.nodes.has(currentRef), objects.has(currentRef));
            const obj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
            cells.set(currentRef, {
                ref: currentRef,
                obj,
                status: 'ok',
                version: transform.version,
                state: { ...ctx.stateCtx.defaultCellState, ...transform.cellState }
            });
            return { action: 'created', obj };
        } else {
            // console.log('updating...', transform.transformer.id);
            const current = cells.get(currentRef)!;
            const oldParams = oldTree.getValue(currentRef)!.params;

            const updateKind = current.status === 'ok' || current.ref === ctx.tree.rootRef
                ? await updateObject(ctx, currentRef, transform.transformer, parent, current.obj!, oldParams, transform.params)
                : Transformer.UpdateResult.Recreate;

            switch (updateKind) {
                case Transformer.UpdateResult.Recreate: {
                    const obj = await createObject(ctx, currentRef, transform.transformer, parent, transform.params);
                    cells.set(currentRef, {
                        ref: currentRef,
                        obj,
                        status: 'ok',
                        version: transform.version,
                        state: { ...ctx.stateCtx.defaultCellState, ...current.state, ...transform.cellState }
                    });
                    return { action: 'replaced', oldObj: current.obj!, newObj: obj };
                }
                case Transformer.UpdateResult.Updated:
                    current.version = transform.version;
                    current.state = { ...ctx.stateCtx.defaultCellState, ...current.state, ...transform.cellState };
                    return { action: 'updated', obj: current.obj };
                default:
                    // TODO check if props need to be updated
                    return { action: 'none' };
            }
        }
    }

    function runTask<T>(t: T | Task<T>, ctx: RuntimeContext) {
        if (typeof (t as any).run === 'function') return (t as Task<T>).runInContext(ctx);
        return t as T;
    }

    function createObject(ctx: UpdateContext, ref: Ref, transformer: Transformer, a: StateObject, params: any) {
        const cache = { };
        ctx.transformCache.set(ref, cache);
        return runTask(transformer.definition.apply({ a, params, cache }, ctx.stateCtx.globalContext), ctx.taskCtx);
    }

    async function updateObject(ctx: UpdateContext, ref: Ref, transformer: Transformer, a: StateObject, b: StateObject, oldParams: any, newParams: any) {
        if (!transformer.definition.update) {
            return Transformer.UpdateResult.Recreate;
        }
        let cache = ctx.transformCache.get(ref);
        if (!cache) {
            cache = { };
            ctx.transformCache.set(ref, cache);
        }
        return runTask(transformer.definition.update({ a, oldParams, b, newParams, cache }, ctx.stateCtx.globalContext), ctx.taskCtx);
    }