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
    readonly nodes: ImmutableTree.Nodes,
    readonly root: ImmutableTree.Node,
    get(ref: ImmutableTree.Ref): Transform | undefined,
    //getChildren(ref: ImmutableTree.Ref): OrderedSet<ImmutableTree.Node>
    asTransient(): TransientTree
}

namespace ImmutableTree {
    export type Ref = Transform.Ref
    export interface Node extends Readonly<TransientTree.Node> { }
    export interface Nodes extends ImmutableMap<ImmutableTree.Ref, Node> { }

    class Impl implements ImmutableTree {
        readonly version: number;
        readonly nodes: ImmutableTree.Nodes;

        get root() { return this.nodes.get(Transform.RootRef)! }

        get(ref: Ref) {
            const n = this.nodes.get(ref);
            return n ? n.value : void 0;
        }

        asTransient(): TransientTree {
            return new TransientTree(this);
        }

        constructor(nodes: ImmutableTree.Nodes, version: number) {
            this.nodes = nodes;
            this.version = version;
        }
    }

    /**
     * Create an instance of an immutable tree.
     */
    export function create<T>(root: Transform): ImmutableTree {
        const r: Node = { value: root, version: 0, parent: root.ref, children: OrderedSet() };
        return new Impl(ImmutableMap([[root.ref, r]]), 0);
    }

    export function fromNodes<T>(nodes: Nodes, version: number): ImmutableTree {
        return new Impl(nodes, version);
    }

    type VisitorCtx = { nodes: Nodes, state: any, f: (node: Node, nodes: Nodes, state: any) => boolean | undefined | void };

    function _postOrderFunc(this: VisitorCtx, c: ImmutableTree.Ref | undefined) { _doPostOrder(this, this.nodes.get(c!)!); }
    function _doPostOrder(ctx: VisitorCtx, root: Node) {
        if (root.children.size) {
            root.children.forEach(_postOrderFunc, ctx);
        }
        ctx.f(root, ctx.nodes, ctx.state);
    }

    /**
     * Visit all nodes in a subtree in "post order", meaning leafs get visited first.
     */
    export function doPostOrder<S>(tree: ImmutableTree, root: Node, state: S, f: (node: Node, nodes: Nodes, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { nodes: tree.nodes, state, f };
        _doPostOrder(ctx, root);
        return ctx.state;
    }

    function _preOrderFunc(this: VisitorCtx, c: ImmutableTree.Ref | undefined) { _doPreOrder(this, this.nodes.get(c!)!); }
    function _doPreOrder(ctx: VisitorCtx, root: Node) {
        const ret = ctx.f(root, ctx.nodes, ctx.state);
        if (typeof ret === 'boolean' && !ret) return;
        if (root.children.size) {
            root.children.forEach(_preOrderFunc, ctx);
        }
    }

    /**
     * Visit all nodes in a subtree in "pre order", meaning leafs get visited last.
     * If the visitor function returns false, the visiting for that branch is interrupted.
     */
    export function doPreOrder<S>(tree: ImmutableTree, root: Node, state: S, f: (node: Node, nodes: Nodes, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { nodes: tree.nodes, state, f };
        _doPreOrder(ctx, root);
        return ctx.state;
    }

    function _subtree(n: Node, nodes: Nodes, subtree: Node[]) { subtree.push(n); }
    /**
     * Get all nodes in a subtree, leafs come first.
     */
    export function subtreePostOrder<T>(tree: ImmutableTree, root: Node) {
        return doPostOrder<Node[]>(tree, root, [], _subtree);
    }


    function _visitChildToJson(this: Ref[], ref: Ref) { this.push(ref); }
    interface ToJsonCtx { nodes: [Transform.Serialized, any, any[]][] }
    function _visitNodeToJson(this: ToJsonCtx, node: Node) {
        const children: Ref[] = [];
        node.children.forEach(_visitChildToJson as any, children);
        this.nodes.push([Transform.toJSON(node.value), node.parent, children]);
    }

    export interface Serialized {
        nodes: [any /** value */, number /** parent index */, number[] /** children indices */][]
    }

    export function toJSON<T>(tree: ImmutableTree): Serialized {
        const ctx: ToJsonCtx = { nodes: [] };

        tree.nodes.forEach(_visitNodeToJson as any, ctx);

        const map = new Map<string, number>();
        let i = 0;
        for (const n of ctx.nodes) map.set(n[0].ref, i++);

        for (const n of ctx.nodes) {
            n[1] = map.get(n[1]);
            const children = n[2];
            for (i = 0; i < children.length; i++) {
                children[i] = map.get(children[i]);
            }
        }
        return {
            nodes: ctx.nodes
        };
    }

    export function fromJSON<T>(data: Serialized): ImmutableTree {
        const nodes = ImmutableMap<ImmutableTree.Ref, Node>().asMutable();

        const values = data.nodes.map(n => Transform.fromJSON(n[0]));
        let i = 0;
        for (const value of values) {
            const node = data.nodes[i++];
            const ref = value.ref;
            nodes.set(ref, {
                value,
                version: 0,
                parent: values[node[1]].ref,
                children: OrderedSet(node[2].map(c => values[c].ref))
            });
        }
        return new Impl(nodes.asImmutable(), 0);
    }
}

class TransientTree implements ImmutableTree {
    nodes = this.tree.nodes.asMutable();
    version: number = this.tree.version + 1;
    private mutations: Map<Transform.Ref, TransientTree.Node> = new Map();

    get root() { return this.nodes.get(Transform.RootRef)! }

    get(ref: Transform.Ref) {
        const n = this.nodes.get(ref);
        return n ? n.value : void 0;
    }

    asTransient() {
        return this.asImmutable().asTransient();
    }

    mutate(ref: ImmutableTree.Ref): TransientTree.Node {
        return mutateNode(this.nodes, this.mutations, ref);
    }

    add(parentRef: ImmutableTree.Ref, value: Transform) {
        const ref = value.ref;
        ensureNotPresent(this.nodes, ref);
        const parent = this.mutate(parentRef);
        const node: TransientTree.Node = { version: 0, value, parent: parentRef, children: OrderedSet<string>().asMutable() };
        this.mutations.set(ref, node);
        parent.children.add(ref);
        this.nodes.set(ref, node);
        return node;
    }

    setValue(ref: ImmutableTree.Ref, value: Transform): TransientTree.Node {
        checkSetRef(ref, value.ref);
        const node = this.mutate(ref);
        node.value = value;
        return node;
    }

    remove(ref: ImmutableTree.Ref): TransientTree.Node[] {
        if (ref === Transform.RootRef) {
            return this.removeChildren(ref);
        }

        const { nodes, mutations } = this;
        const node = nodes.get(ref);
        if (!node) return [];
        this.mutate(node.parent).children.delete(ref);

        const st = ImmutableTree.subtreePostOrder(this, node);
        for (const n of st) {
            nodes.delete(n.value.ref);
            mutations.delete(n.value.ref);
        }

        return st;
    }

    removeChildren(ref: ImmutableTree.Ref): TransientTree.Node[] {
        const { nodes, mutations } = this;
        let node = nodes.get(ref);
        if (!node || !node.children.size) return [];
        node = this.mutate(ref);
        const st = ImmutableTree.subtreePostOrder(this, node);
        // remove the last node which is the parent
        st.pop();
        node.children.clear();
        for (const n of st) {
            nodes.delete(n.value.ref);
            mutations.delete(n.value.ref);
        }
        return st;
    }

    asImmutable() {
        if (this.mutations.size === 0) return this.tree;

        this.mutations.forEach(m => (m as TransientTree.Node).children = m.children.asImmutable());
        return ImmutableTree.fromNodes(this.nodes.asMutable(), this.version);
    }

    constructor(private tree: ImmutableTree) {

    }
}

namespace TransientTree {
    export interface Node { value: Transform, version: number, parent: ImmutableTree.Ref, children: OrderedSet<ImmutableTree.Ref> }
}


function checkSetRef(oldRef: ImmutableTree.Ref, newRef: ImmutableTree.Ref) {
    if (oldRef !== newRef) {
        throw new Error(`Cannot setValue of node '${oldRef}' because the new value has a different ref '${newRef}'.`);
    }
}

function ensureNotPresent(nodes: ImmutableTree.Nodes, ref: ImmutableTree.Ref) {
    if (nodes.has(ref)) {
        throw new Error(`Cannot add node '${ref}' because a different node with this ref already present in the tree.`);
    }
}

function ensurePresent(nodes: ImmutableTree.Nodes, ref: ImmutableTree.Ref) {
    if (!nodes.has(ref)) {
        throw new Error(`Node '${ref}' is not present in the tree.`);
    }
}

function mutateNode(nodes: ImmutableTree.Nodes, mutations: Map<ImmutableTree.Ref, TransientTree.Node>, ref: ImmutableTree.Ref): TransientTree.Node {
    ensurePresent(nodes, ref);
    if (mutations.has(ref)) {
        return mutations.get(ref)!;
    }
    const node = nodes.get(ref)!;
    const newNode: TransientTree.Node = { value: node.value, version: node.version + 1, parent: node.parent, children: node.children.asMutable() };
    mutations.set(ref, newNode);
    nodes.set(ref, newNode);
    return newNode;
}