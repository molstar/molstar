/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';
import { StateTransform } from '../transform';
import { TransientTree } from './transient';

export { StateTree };

/**
 * An immutable tree where each node requires a unique reference.
 * Represented as an immutable map.
 */
interface StateTree {
    readonly root: StateTransform,
    readonly transforms: StateTree.Transforms,
    readonly children: StateTree.Children,
    /** Refs to all nodes that depend on the given key */
    readonly dependencies: StateTree.Dependencies,

    asTransient(): TransientTree
}

namespace StateTree {
    type Ref = StateTransform.Ref

    export interface ChildSet {
        readonly size: number,
        readonly values: OrderedSet<Ref>['values'],
        has(ref: Ref): boolean,
        readonly forEach: OrderedSet<Ref>['forEach'],
        readonly map: OrderedSet<Ref>['map'],
        toArray(): Ref[],
        first(): Ref
    }

    interface _Map<T> {
        readonly size: number,
        has(ref: Ref): boolean,
        get(ref: Ref): T
    }

    export interface Transforms extends _Map<StateTransform> {}
    export interface Children extends _Map<ChildSet> { }
    export interface Dependencies extends _Map<ChildSet> { }

    class Impl implements StateTree {
        get root() { return this.transforms.get(StateTransform.RootRef)!; }

        asTransient(): TransientTree {
            return new TransientTree(this);
        }

        constructor(public transforms: Transforms, public children: Children, public dependencies: Dependencies) {
        }
    }

    /**
     * Create an instance of an immutable tree.
     */
    export function createEmpty(customRoot?: StateTransform): StateTree {
        const root = customRoot || StateTransform.createRoot();
        return create(ImmutableMap([[root.ref, root]]), ImmutableMap([[root.ref, OrderedSet()]]), ImmutableMap());
    }

    export function create(nodes: Transforms, children: Children, dependencies: Dependencies): StateTree {
        return new Impl(nodes, children, dependencies);
    }

    type VisitorCtx = { tree: StateTree, state: any, f: (node: StateTransform, tree: StateTree, state: any) => boolean | undefined | void };

    function _postOrderFunc(this: VisitorCtx, c: Ref | undefined) { _doPostOrder(this, this.tree.transforms.get(c!)!); }
    function _doPostOrder(ctx: VisitorCtx, root: StateTransform) {
        const children = ctx.tree.children.get(root.ref);
        if (children && children.size) {
            children.forEach(_postOrderFunc, ctx);
        }
        ctx.f(root, ctx.tree, ctx.state);
    }

    /**
     * Visit all nodes in a subtree in "post order", meaning leafs get visited first.
     */
    export function doPostOrder<S>(tree: StateTree, root: StateTransform, state: S, f: (node: StateTransform, tree: StateTree, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { tree, state, f };
        _doPostOrder(ctx, root);
        return ctx.state;
    }

    function _preOrderFunc(this: VisitorCtx, c: Ref | undefined) { _doPreOrder(this, this.tree.transforms.get(c!)!); }
    function _doPreOrder(ctx: VisitorCtx, root: StateTransform) {
        const ret = ctx.f(root, ctx.tree, ctx.state);
        if (typeof ret === 'boolean' && !ret) return;
        const children = ctx.tree.children.get(root.ref);
        if (children && children.size) {
            children.forEach(_preOrderFunc, ctx);
        }
    }

    /**
     * Visit all nodes in a subtree in "pre order", meaning leafs get visited last.
     * If the visitor function returns false, the visiting for that branch is interrupted.
     */
    export function doPreOrder<S>(tree: StateTree, root: StateTransform, state: S, f: (node: StateTransform, tree: StateTree, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { tree, state, f };
        _doPreOrder(ctx, root);
        return ctx.state;
    }

    function _subtree(n: StateTransform, _: any, subtree: StateTransform[]) { subtree.push(n); }
    /**
     * Get all nodes in a subtree, leafs come first.
     */
    export function subtreePostOrder(tree: StateTree, root: StateTransform) {
        return doPostOrder<StateTransform[]>(tree, root, [], _subtree);
    }

    function _visitNodeToJson(node: StateTransform, tree: StateTree, ctx: StateTransform.Serialized[]) {
        // const children: Ref[] = [];
        // tree.children.get(node.ref).forEach(_visitChildToJson as any, children);
        ctx.push(StateTransform.toJSON(node));
    }

    export interface Serialized {
        /** Transforms serialized in pre-order */
        transforms: StateTransform.Serialized[]
    }

    export function toJSON(tree: StateTree): Serialized {
        const transforms: StateTransform.Serialized[] = [];
        doPreOrder(tree, tree.root, transforms, _visitNodeToJson);
        return { transforms };
    }

    export function fromJSON(data: Serialized): StateTree {
        const nodes = ImmutableMap<Ref, StateTransform>().asMutable();
        const children = ImmutableMap<Ref, OrderedSet<Ref>>().asMutable();
        const dependencies = ImmutableMap<Ref, OrderedSet<Ref>>().asMutable();

        for (const t of data.transforms) {
            const transform = StateTransform.fromJSON(t);
            nodes.set(transform.ref, transform);

            if (!children.has(transform.ref)) {
                children.set(transform.ref, OrderedSet<Ref>().asMutable());
            }

            if (transform.ref !== transform.parent) children.get(transform.parent).add(transform.ref);
        }

        const dependent = new Set<Ref>();
        for (const t of data.transforms) {
            const ref = t.ref;
            children.set(ref, children.get(ref).asImmutable());

            if (!t.dependsOn) continue;

            for (const d of t.dependsOn) {
                dependent.add(d);
                if (!dependencies.has(d)) {
                    dependencies.set(d, OrderedSet<Ref>([ref]).asMutable());
                } else {
                    dependencies.get(d).add(ref);
                }
            }
        }

        dependent.forEach(d => {
            dependencies.set(d, dependencies.get(d).asImmutable());
        });

        return create(nodes.asImmutable(), children.asImmutable(), dependencies.asImmutable());
    }

    export function dump(tree: StateTree) {
        console.log({
            tr: (tree.transforms as ImmutableMap<any, any>).keySeq().toArray(),
            tr1: (tree.transforms as ImmutableMap<any, any>).valueSeq().toArray().map(t => t.ref),
            ch: (tree.children as ImmutableMap<any, any>).keySeq().toArray()
        });
    }

    function _subtreeHasRef(tree: StateTree, root: StateTransform.Ref, ref: StateTransform.Ref) {
        if (root === ref) return true;
        const children = tree.children.get(root);
        const it = children.values();
        while (true) {
            const next = it.next();
            if (next.done) return false;
            if (_subtreeHasRef(tree, next.value, ref)) return true;
        }
    }

    /** Check if the subtree with the given root has the provided ref */
    export function subtreeHasRef(tree: StateTree, root: StateTransform.Ref, ref: StateTransform.Ref): boolean {
        if (!tree.transforms.has(root) || !tree.transforms.has(ref)) return false;
        return _subtreeHasRef(tree, root, ref);
    }

    export function getDecoratorRoot(tree: StateTree, ref: StateTransform.Ref): StateTransform.Ref {
        const children = tree.children.get(ref);
        if (children.size !== 1) return ref;
        const child = tree.transforms.get(children.first());
        if (child.transformer.definition.isDecorator) return getDecoratorRoot(tree, child.ref);
        return ref;
    }
}