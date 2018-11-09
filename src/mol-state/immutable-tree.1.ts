/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';
import { Transform } from './transform';

export { ImmutableTree, TransientTree }

/**
 * An immutable tree where each node requires a unique reference.
 * Represented as an immutable map.
 */
interface ImmutableTree {
    readonly version: number,
    readonly root: Transform,
    readonly nodes: ImmutableTree.Nodes,
    readonly children: ImmutableTree.Children,
    asTransient(): TransientTree
}

namespace ImmutableTree {
    type Ref = Transform.Ref

    export interface ChildSet {
        readonly size: number,
        readonly values: OrderedSet<Ref>['values'],
        has(ref: Ref): boolean,
        readonly forEach: OrderedSet<Ref>['forEach']
    }

    export type Node = Transform
    export type Nodes = ImmutableMap<Ref, Transform>
    export type Children = ImmutableMap<Ref, ChildSet>

    class Impl implements ImmutableTree {
        get root() { return this.nodes.get(Transform.RootRef)! }

        asTransient(): TransientTree {
            return new TransientTree(this);
        }

        constructor(public nodes: ImmutableTree.Nodes, public children: Children, public version: number) {
        }
    }

    /**
     * Create an instance of an immutable tree.
     */
    export function create<T>(root: Transform): ImmutableTree {
        return new Impl(ImmutableMap([[root.ref, root]]), ImmutableMap([[root.ref, OrderedSet()]]), 0);
    }

    export function construct<T>(nodes: Nodes, children: Children, version: number): ImmutableTree {
        return new Impl(nodes, children, version);
    }

    type VisitorCtx = { tree: ImmutableTree, state: any, f: (node: Node, tree: ImmutableTree, state: any) => boolean | undefined | void };

    function _postOrderFunc(this: VisitorCtx, c: Ref | undefined) { _doPostOrder(this, this.tree.nodes.get(c!)); }
    function _doPostOrder(ctx: VisitorCtx, root: Node) {
        const children = ctx.tree.children.get(root.ref);
        if (checkSetRef && children.size) {
            children.forEach(_postOrderFunc, ctx);
        }
        ctx.f(root, ctx.tree, ctx.state);
    }

    /**
     * Visit all nodes in a subtree in "post order", meaning leafs get visited first.
     */
    export function doPostOrder<S>(tree: ImmutableTree, root: Node, state: S, f: (node: Node, tree: ImmutableTree, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { tree, state, f };
        _doPostOrder(ctx, root);
        return ctx.state;
    }

    function _preOrderFunc(this: VisitorCtx, c: Ref | undefined) { _doPreOrder(this, this.tree.nodes.get(c!)); }
    function _doPreOrder(ctx: VisitorCtx, root: Node) {
        const ret = ctx.f(root, ctx.tree, ctx.state);
        if (typeof ret === 'boolean' && !ret) return;
        const children = ctx.tree.children.get(root.ref);
        if (checkSetRef && children.size) {
            children.forEach(_preOrderFunc, ctx);
        }
    }

    /**
     * Visit all nodes in a subtree in "pre order", meaning leafs get visited last.
     * If the visitor function returns false, the visiting for that branch is interrupted.
     */
    export function doPreOrder<S>(tree: ImmutableTree, root: Node, state: S, f: (node: Node, tree: ImmutableTree, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { tree, state, f };
        _doPreOrder(ctx, root);
        return ctx.state;
    }

    function _subtree(n: Node, _: any, subtree: Node[]) { subtree.push(n); }
    /**
     * Get all nodes in a subtree, leafs come first.
     */
    export function subtreePostOrder<T>(tree: ImmutableTree, root: Node) {
        return doPostOrder<Node[]>(tree, root, [], _subtree);
    }


    function _visitChildToJson(this: Ref[], ref: Ref) { this.push(ref); }
    interface ToJsonCtx { tree: ImmutableTree, nodes: [Transform.Serialized, any[]][] }
    function _visitNodeToJson(this: ToJsonCtx, node: Node) {
        const children: Ref[] = [];
        this.tree.children.get(node.ref).forEach(_visitChildToJson as any, children);
        this.nodes.push([Transform.toJSON(node), children]);
    }

    export interface Serialized {
        nodes: [any /** value */, number[] /** children indices */][]
    }

    export function toJSON<T>(tree: ImmutableTree): Serialized {
        const ctx: ToJsonCtx = { tree, nodes: [] };

        tree.nodes.forEach(_visitNodeToJson as any, ctx);

        const map = new Map<string, number>();
        let i = 0;
        for (const n of ctx.nodes) map.set(n[0].ref, i++);

        for (const n of ctx.nodes) {
            const children = n[1];
            for (i = 0; i < children.length; i++) {
                children[i] = map.get(children[i]);
            }
        }
        return {
            nodes: ctx.nodes
        };
    }

    export function fromJSON<T>(data: Serialized): ImmutableTree {
        const nodes = ImmutableMap<Ref, Node>().asMutable();
        const children = ImmutableMap<Ref, OrderedSet<Ref>>().asMutable();

        const values = data.nodes.map(n => Transform.fromJSON(n[0]));
        let i = 0;
        for (const value of values) {
            const node = data.nodes[i++];
            const ref = value.ref;
            nodes.set(ref, value);
            children.set(ref, OrderedSet(node[1].map(c => values[c].ref)));
        }
        return new Impl(nodes.asImmutable(), children.asImmutable(), 0);
    }
}

class TransientTree implements ImmutableTree {
    nodes = this.tree.nodes.asMutable();
    children = this.tree.children.asMutable();

    version: number = this.tree.version + 1;

    private mutations: Map<Transform.Ref, OrderedSet<Transform.Ref>> = new Map();

    get root() { return this.nodes.get(Transform.RootRef)! }

    asTransient() {
        return this.asImmutable().asTransient();
    }

    private addChild(parent: Transform.Ref, child: Transform.Ref) {
        if (this.mutations.has(parent)) {
            this.mutations.get(parent)!.add(child);
        } else {
            const set = (this.children.get(parent) as OrderedSet<Transform.Ref>).asMutable();
            set.add(child);
            this.children.set(parent, set);
            this.mutations.set(parent, set);
        }
    }

    private removeChild(parent: Transform.Ref, child: Transform.Ref) {
        if (this.mutations.has(parent)) {
            this.mutations.get(parent)!.remove(child);
        } else {
            const set = (this.children.get(parent) as OrderedSet<Transform.Ref>).asMutable();
            set.remove(child);
            this.children.set(parent, set);
            this.mutations.set(parent, set);
        }
    }

    private clearRoot() {
        const parent = Transform.RootRef;
        if (this.children.get(parent).size === 0) return;
        const set = OrderedSet<Transform.Ref>();
        this.children.set(parent, set);
        this.mutations.set(parent, set);
    }

    add(transform: Transform) {
        const ref = transform.ref;

        const children = this.children.get(transform.parent);
        if (!children) parentNotPresent(transform.parent);

        if (!children.has(transform.ref)) {
            this.addChild(transform.parent, transform.ref);
        }

        if (!this.children.has(transform.ref)) {
            this.children.set(transform.ref, OrderedSet());
        }

        this.nodes.set(ref, transform);
        return transform;
    }

    set(transform: Transform) {
        ensurePresent(this.nodes, transform.ref);
        this.nodes.set(transform.ref, transform);
    }

    remove(ref: Transform.Ref): Transform[] {
        const { nodes, mutations, children } = this;
        const node = nodes.get(ref);
        if (!node) return [];

        const st = ImmutableTree.subtreePostOrder(this, node);
        if (ref === Transform.RootRef) {
            st.pop();
            this.clearRoot();
        } else {
            this.removeChild(node.parent, node.ref);
        }

        for (const n of st) {
            nodes.delete(n.ref);
            children.delete(n.ref);
            mutations.delete(n.ref);
        }

        return st;
    }

    // removeChildren(ref: ImmutableTree.Ref): TransientTree.Node[] {
    //     const { nodes, mutations } = this;
    //     let node = nodes.get(ref);
    //     if (!node || !node.children.size) return [];
    //     node = this.mutate(ref);
    //     const st = ImmutableTree.subtreePostOrder(this, node);
    //     // remove the last node which is the parent
    //     st.pop();
    //     node.children.clear();
    //     for (const n of st) {
    //         nodes.delete(n.value.ref);
    //         mutations.delete(n.value.ref);
    //     }
    //     return st;
    // }

    asImmutable() {
        if (this.mutations.size === 0) return this.tree;

        this.mutations.forEach(fixChildMutations, this.children);
        return ImmutableTree.construct(this.nodes.asImmutable(), this.children.asImmutable(), this.version);
    }

    constructor(private tree: ImmutableTree) {

    }
}

function fixChildMutations(this: ImmutableTree.Children, m: OrderedSet<Transform.Ref>, k: Transform.Ref) { this.set(k, m.asImmutable()); }

function checkSetRef(oldRef: Transform.Ref, newRef: Transform.Ref) {
    if (oldRef !== newRef) {
        throw new Error(`Cannot setValue of node '${oldRef}' because the new value has a different ref '${newRef}'.`);
    }
}

function parentNotPresent(ref: Transform.Ref) {
    throw new Error(`Parent '${ref}' must be present in the tree.`);
}

function ensurePresent(nodes: ImmutableTree.Nodes, ref: Transform.Ref) {
    if (!nodes.has(ref)) {
        throw new Error(`Node '${ref}' is not present in the tree.`);
    }
}