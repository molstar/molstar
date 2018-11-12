/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';
import { Transform } from '../transform';
import { TransientTree } from './transient';
import { StateTreeBuilder } from './builder';

export { StateTree }

/**
 * An immutable tree where each node requires a unique reference.
 * Represented as an immutable map.
 */
interface StateTree {
    readonly root: Transform,
    readonly transforms: StateTree.Nodes,
    readonly children: StateTree.Children,
    asTransient(): TransientTree,
    build(): StateTreeBuilder.Root
}

namespace StateTree {
    type Ref = Transform.Ref

    export interface ChildSet {
        readonly size: number,
        readonly values: OrderedSet<Ref>['values'],
        has(ref: Ref): boolean,
        readonly forEach: OrderedSet<Ref>['forEach'],
        readonly map: OrderedSet<Ref>['map']
    }

    interface _Map<T> {
        readonly size: number,
        has(ref: Ref): boolean,
        get(ref: Ref): T
    }

    export interface Nodes extends _Map<Transform> {}
    export interface Children extends _Map<ChildSet> { }

    class Impl implements StateTree {
        get root() { return this.transforms.get(Transform.RootRef)! }

        asTransient(): TransientTree {
            return new TransientTree(this);
        }

        build(): StateTreeBuilder.Root {
            return new StateTreeBuilder.Root(this);
        }

        constructor(public transforms: StateTree.Nodes, public children: Children) {
        }
    }

    /**
     * Create an instance of an immutable tree.
     */
    export function createEmpty(): StateTree {
        const root = Transform.createRoot();
        return create(ImmutableMap([[root.ref, root]]), ImmutableMap([[root.ref, OrderedSet()]]));
    }

    export function create(nodes: Nodes, children: Children): StateTree {
        return new Impl(nodes, children);
    }

    type VisitorCtx = { tree: StateTree, state: any, f: (node: Transform, tree: StateTree, state: any) => boolean | undefined | void };

    function _postOrderFunc(this: VisitorCtx, c: Ref | undefined) { _doPostOrder(this, this.tree.transforms.get(c!)!); }
    function _doPostOrder(ctx: VisitorCtx, root: Transform) {
        const children = ctx.tree.children.get(root.ref);
        if (children && children.size) {
            children.forEach(_postOrderFunc, ctx);
        }
        ctx.f(root, ctx.tree, ctx.state);
    }

    /**
     * Visit all nodes in a subtree in "post order", meaning leafs get visited first.
     */
    export function doPostOrder<S>(tree: StateTree, root: Transform, state: S, f: (node: Transform, tree: StateTree, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { tree, state, f };
        _doPostOrder(ctx, root);
        return ctx.state;
    }

    function _preOrderFunc(this: VisitorCtx, c: Ref | undefined) { _doPreOrder(this, this.tree.transforms.get(c!)!); }
    function _doPreOrder(ctx: VisitorCtx, root: Transform) {
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
    export function doPreOrder<S>(tree: StateTree, root: Transform, state: S, f: (node: Transform, tree: StateTree, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { tree, state, f };
        _doPreOrder(ctx, root);
        return ctx.state;
    }

    function _subtree(n: Transform, _: any, subtree: Transform[]) { subtree.push(n); }
    /**
     * Get all nodes in a subtree, leafs come first.
     */
    export function subtreePostOrder<T>(tree: StateTree, root: Transform) {
        return doPostOrder<Transform[]>(tree, root, [], _subtree);
    }


    // function _visitChildToJson(this: Ref[], ref: Ref) { this.push(ref); }
    // interface ToJsonCtx { nodes: Transform.Serialized[] }
    function _visitNodeToJson(node: Transform, tree: StateTree, ctx: Transform.Serialized[]) {
        // const children: Ref[] = [];
        // tree.children.get(node.ref).forEach(_visitChildToJson as any, children);
        ctx.push(Transform.toJSON(node));
    }

    export interface Serialized {
        /** Tree nodes serialized in pre-order */
        nodes: Transform.Serialized[]
    }

    export function toJSON<T>(tree: StateTree): Serialized {
        const nodes: Transform.Serialized[] = [];
        doPreOrder(tree, tree.root, nodes, _visitNodeToJson);
        return { nodes };
    }

    export function fromJSON<T>(data: Serialized): StateTree {
        const nodes = ImmutableMap<Ref, Transform>().asMutable();
        const children = ImmutableMap<Ref, OrderedSet<Ref>>().asMutable();

        for (const s of data.nodes) {
            const node = Transform.fromJSON(s);
            nodes.set(node.ref, node);

            if (!children.has(node.ref)) {
                children.set(node.ref, OrderedSet<Ref>().asMutable());
            }

            if (node.ref !== node.parent) children.get(node.parent).add(node.ref);
        }

        for (const s of data.nodes) {
            children.set(s.ref, children.get(s.ref).asImmutable());
        }

        return create(nodes.asImmutable(), children.asImmutable());
    }
}