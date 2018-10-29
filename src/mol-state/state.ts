/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { StateTree } from './tree';
import { Transform } from './tree/transform';
import { ImmutableTree } from './util/immutable-tree';
import { Transformer } from './transformer';
import { StateContext } from './context';
import { UUID } from 'mol-util';
import { RuntimeContext, Task } from 'mol-task';

export interface State {
    tree: StateTree,
    objects: State.Objects,
    context: StateContext
}

export namespace State {
    export type Ref = Transform.Ref
    export type Objects = Map<Ref, StateObject.Node>

    export function create(params?: { globalContext?: unknown, defaultObjectProps: unknown }) {
        const tree = StateTree.create();
        const objects: Objects = new Map();
        const root = tree.getValue(tree.rootRef)!;
        const defaultObjectProps = (params && params.defaultObjectProps) || { }

        objects.set(tree.rootRef, { obj: void 0 as any, state: StateObject.StateType.Ok, version: root.version, props: { ...defaultObjectProps } });

        return {
            tree,
            objects,
            context: StateContext.create({
                globalContext: params && params.globalContext,
                defaultObjectProps
            })
        };
    }

    export function update(state: State, tree: StateTree): Task<State> {
        return Task.create('Update Tree', taskCtx => {
            const ctx: UpdateContext = {
                stateCtx: state.context,
                taskCtx,
                oldTree: state.tree,
                tree: tree,
                objects: state.objects
            };
            return _update(ctx);
        })
    }

    async function _update(ctx: UpdateContext): Promise<State> {
        const roots = findUpdateRoots(ctx.objects, ctx.tree);
        const deletes = findDeletes(ctx);
        for (const d of deletes) {
            ctx.objects.delete(d);
            ctx.stateCtx.events.object.removed.next({ ref: d });
        }

        initObjectState(ctx, roots);

        for (const root of roots) {
            await updateSubtree(ctx, root);
        }

        return {
            tree: ctx.tree,
            objects: ctx.objects,
            context: ctx.stateCtx
        };
    }

    interface UpdateContext {
        stateCtx: StateContext,
        taskCtx: RuntimeContext,
        oldTree: StateTree,
        tree: StateTree,
        objects: Objects
    }

    function findUpdateRoots(objects: Objects, tree: StateTree) {
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

    function findDeletes(ctx: UpdateContext): Ref[] {
        // TODO: do this in some sort of "tree order"?
        const deletes: Ref[] = [];
        const keys = ctx.objects.keys();
        while (true) {
            const key = keys.next();
            if (key.done) break;
            if (!ctx.tree.nodes.has(key.value)) deletes.push(key.value);
        }
        return deletes;
    }

    function setObjectState(ctx: UpdateContext, ref: Ref, state: StateObject.StateType, errorText?: string) {
        let changed = false;
        if (ctx.objects.has(ref)) {
            const obj = ctx.objects.get(ref)!;
            changed = obj.state !== state;
            obj.state = state;
            obj.errorText = errorText;
        } else {
            const obj = { state, version: UUID.create(), errorText, props: { ...ctx.stateCtx.defaultObjectProps } };
            ctx.objects.set(ref, obj);
            changed = true;
        }
        if (changed) ctx.stateCtx.events.object.stateChanged.next({ ref });
    }

    function _initVisitor(t: ImmutableTree.Node<Transform>, _: any, ctx: UpdateContext) {
        setObjectState(ctx, t.ref, StateObject.StateType.Pending);
    }
    /** Return "resolve set" */
    function initObjectState(ctx: UpdateContext, roots: Ref[]) {
        for (const root of roots) {
            ImmutableTree.doPreOrder(ctx.tree, ctx.tree.nodes.get(root), ctx, _initVisitor);
        }
    }

    function doError(ctx: UpdateContext, ref: Ref, errorText: string) {
        setObjectState(ctx, ref, StateObject.StateType.Error, errorText);
        const wrap = ctx.objects.get(ref)!;
        if (wrap.obj) {
            ctx.stateCtx.events.object.removed.next({ ref });
            wrap.obj = void 0;
        }

        const children = ctx.tree.nodes.get(ref)!.children.values();
        while (true) {
            const next = children.next();
            if (next.done) return;
            doError(ctx, next.value, 'Parent node contains error.');
        }
    }

    function findParent(tree: StateTree, objects: Objects, root: Ref, types: { type: StateObject.Type }[]): StateObject {
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
        setObjectState(ctx, root, StateObject.StateType.Processing);

        try {
            const update = await updateNode(ctx, root);
            setObjectState(ctx, root, StateObject.StateType.Ok);
            if (update.action === 'created') {
                ctx.stateCtx.events.object.created.next({ ref: root });
            } else if (update.action === 'updated') {
                ctx.stateCtx.events.object.updated.next({ ref: root });
            } else if (update.action === 'replaced') {
                ctx.stateCtx.events.object.replaced.next({ ref: root, old: update.old });
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
        const { oldTree, tree, objects } = ctx;
        const transform = tree.getValue(currentRef)!;
        const parent = findParent(tree, objects, currentRef, transform.transformer.definition.from);
        // console.log('parent', transform.transformer.id, transform.transformer.definition.from[0].type, parent ? parent.ref : 'undefined')
        if (!oldTree.nodes.has(currentRef) || !objects.has(currentRef)) {
            // console.log('creating...', transform.transformer.id, oldTree.nodes.has(currentRef), objects.has(currentRef));
            const obj = await createObject(ctx, transform.transformer, parent, transform.params);
            obj.ref = currentRef;
            objects.set(currentRef, {
                obj,
                state: StateObject.StateType.Ok,
                version: transform.version,
                props: { ...ctx.stateCtx.defaultObjectProps, ...transform.defaultProps }
            });
            return { action: 'created' };
        } else {
            // console.log('updating...', transform.transformer.id);
            const current = objects.get(currentRef)!;
            const oldParams = oldTree.getValue(currentRef)!.params;
            switch (await updateObject(ctx, transform.transformer, parent, current.obj!, oldParams, transform.params)) {
                case Transformer.UpdateResult.Recreate: {
                    const obj = await createObject(ctx, transform.transformer, parent, transform.params);
                    obj.ref = currentRef;
                    objects.set(currentRef, {
                        obj,
                        state: StateObject.StateType.Ok,
                        version: transform.version,
                        props: { ...ctx.stateCtx.defaultObjectProps, ...current.props, ...transform.defaultProps }
                    });
                    return { action: 'replaced', old: current.obj! };
                }
                case Transformer.UpdateResult.Updated:
                    current.version = transform.version;
                    current.props = { ...ctx.stateCtx.defaultObjectProps, ...current.props, ...transform.defaultProps };
                    return { action: 'updated' };
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

    function createObject(ctx: UpdateContext, transformer: Transformer, a: StateObject, params: any) {
        return runTask(transformer.definition.apply({ a, params }, ctx.stateCtx.globalContext), ctx.taskCtx);
    }

    async function updateObject(ctx: UpdateContext, transformer: Transformer, a: StateObject, b: StateObject, oldParams: any, newParams: any) {
        if (!transformer.definition.update) {
            return Transformer.UpdateResult.Recreate;
        }
        return runTask(transformer.definition.update({ a, oldParams, b, newParams }, ctx.stateCtx.globalContext), ctx.taskCtx);
    }
}
