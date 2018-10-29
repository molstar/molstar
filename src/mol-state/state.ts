/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { TransformTree } from './tree/tree';
import { Transform } from './tree/transform';
import { Map as ImmutableMap } from 'immutable';
// import { StateContext } from './context/context';
import { ImmutableTree } from './util/immutable-tree';
import { Transformer } from './transformer';
import { Task } from 'mol-task';

export interface State<ObjectProps = unknown> {
    definition: State.Definition<ObjectProps>,
    objects: State.Objects
}

export namespace State {
    export type ObjectProps<P> = ImmutableMap<Transform.Ref, P>
    export type Objects = Map<Transform.Ref, StateObject.Wrapped>

    export interface Definition<P = unknown> {
        tree: TransformTree,
        // things like object visibility
        props: ObjectProps<P>
    }

    export function create(): State {
        const tree = TransformTree.create();
        const objects: Objects = new Map();
        const root = tree.getValue(tree.rootRef)!;

        objects.set(tree.rootRef, { obj: void 0 as any, state: StateObject.StateType.Ok, version: root.version });

        return {
            definition: {
                tree,
                props: ImmutableMap()
            },
            objects
        };
    }

    export async function update<P>(state: State<P>, tree: TransformTree, props?: ObjectProps<P>): Promise<State<P>> {
        const roots = findUpdateRoots(state.objects, tree);
        const deletes = findDeletes(state.objects, tree);
        for (const d of deletes) {
            state.objects.delete(d);
        }

        console.log('roots', roots);
        for (const root of roots) {
            await updateSubtree(state.definition.tree, tree, state.objects, root);
        }

        return {
            definition: { tree, props: props || state.definition.props },
            objects: state.objects
        };
    }

    function findUpdateRoots(objects: Objects, tree: TransformTree) {
        console.log(tree);
        const findState = {
            roots: [] as Transform.Ref[],
            objects
        };

        ImmutableTree.doPreOrder(tree, tree.nodes.get(tree.rootRef)!, findState, (n, _, s) => {
            if (!s.objects.has(n.ref)) {
                console.log('missing', n.ref);
                s.roots.push(n.ref);
                return false;
            }
            const o = s.objects.get(n.ref)!;
            if (o.version !== n.value.version) {
                console.log('diff version', n.ref, n.value.version, o.version);
                s.roots.push(n.ref);
                return false;
            }

            return true;
        });

        return findState.roots;
    }

    function findDeletes(objects: Objects, tree: TransformTree): Transform.Ref[] {
        // TODO
        return [];
    }

    function findParent(tree: TransformTree, objects: Objects, root: Transform.Ref, types: { type: StateObject.Type }[]): StateObject {
        let current = tree.nodes.get(root)!;
        console.log('finding', types.map(t => t.type.kind));
        while (true) {
            current = tree.nodes.get(current.parent)!;
            if (current.ref === tree.rootRef) return objects.get(tree.rootRef)!.obj;
            const obj = objects.get(current.ref)!.obj;
            console.log('current', obj.type.kind);
            for (const t of types) if (obj.type === t.type) return objects.get(current.ref)!.obj;
        }
    }

    async function updateSubtree(oldTree: TransformTree, tree: TransformTree, objects: Objects, root: Transform.Ref) {
        await updateNode(oldTree, tree, objects, root);
        const children = tree.nodes.get(root)!.children.values();
        while (true) {
            const next = children.next();
            if (next.done) return;
            await updateSubtree(oldTree, tree, objects, next.value);
        }
    }

    async function updateNode(oldTree: TransformTree, tree: TransformTree, objects: Objects, currentRef: Transform.Ref) {
        const transform = tree.getValue(currentRef)!;
        const parent = findParent(tree, objects, currentRef, transform.transformer.definition.from);
        console.log('parent', parent ? parent.ref : 'undefined')
        if (!oldTree.nodes.has(transform.ref) || !objects.has(transform.ref)) {
            console.log('creating...', transform.transformer.id, oldTree.nodes.has(transform.ref), objects.has(transform.ref));
            const obj = await createObject(transform.transformer, parent, transform.params);
            obj.ref = transform.ref;
            objects.set(currentRef, { obj, state: StateObject.StateType.Ok, version: transform.version });
        } else {
            console.log('updating...', transform.transformer.id);
            const current = objects.get(transform.ref)!.obj;
            const oldParams = oldTree.getValue(transform.ref)!.params;
            switch (await updateObject(transform.transformer, parent, current, oldParams, transform.params)) {
                case Transformer.UpdateResult.Recreate: {
                    const obj = await createObject(transform.transformer, parent, transform.params);
                    obj.ref = transform.ref;
                    objects.set(currentRef, { obj, state: StateObject.StateType.Ok, version: transform.version });
                    break;
                }
                case Transformer.UpdateResult.Updated: {
                    const obj = objects.get(currentRef)!;
                    obj.version = transform.version;
                    break;
                }
            }
        }
    }

    async function runTask<A>(t: A | Task<A>): Promise<A> {
        if ((t as any).run) return await (t as Task<A>).run();
        return t as A;
    }

    function createObject(transformer: Transformer, a: StateObject, params: any) {
        return runTask(transformer.definition.apply({ a, params }));
    }

    async function updateObject(transformer: Transformer, a: StateObject, b: StateObject, oldParams: any, newParams: any) {
        if (!transformer.definition.update) {
            return Transformer.UpdateResult.Recreate;
        }
        return runTask(transformer.definition.update({ a, oldParams, b, newParams }));
    }
}
