/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';

/**
 * An immutable tree where each node requires a unique reference.
 * Represented as an immutable map.
 */
export interface ImmutableTree<T> {
    readonly rootRef: string,
    readonly version: number,
    readonly nodes: ImmutableTree.Nodes<T>,
    getRef(e: T): string
}

export namespace ImmutableTree {
    export interface MutableNode<T> { ref: string, value: T, version: number, parent: string, children: OrderedSet<string> }
    export interface Node<T> extends Readonly<MutableNode<T>> { }
    export interface Nodes<T> extends ImmutableMap<string, Node<T>> { }

    class Impl<T> implements ImmutableTree<T> {
        readonly rootRef: string;
        readonly version: number;
        readonly nodes: ImmutableTree.Nodes<T>;
        readonly getRef: (e: T) => string;

        constructor(rootRef: string, nodes: ImmutableTree.Nodes<T>, getRef: (e: T) => string, version: number) {
            this.rootRef = rootRef;
            this.nodes = nodes;
            this.getRef = getRef;
            this.version = version;
        }
    }

    /**
     * Create an instance of an immutable tree.
     */
    export function create<T>(root: T, getRef: (t: T) => string): ImmutableTree<T> {
        const ref = getRef(root);
        const r: Node<T> = { ref, value: root, version: 0, parent: ref, children: OrderedSet() };
        return new Impl(ref, ImmutableMap([[ref, r]]), getRef, 0);
    }

    export function asTransient<T>(tree: ImmutableTree<T>) {
        return new Transient(tree);
    }

    type N = Node<any>
    type Ns = Nodes<any>

    type VisitorCtx = { nodes: Ns, state: any, f: (node: N, nodes: Ns, state: any) => boolean | undefined | void };

    function _postOrderFunc(this: VisitorCtx, c: string | undefined) { _doPostOrder(this, this.nodes.get(c!)!); }
    function _doPostOrder<T, S>(ctx: VisitorCtx, root: N) {
        if (root.children.size) {
            root.children.forEach(_postOrderFunc, ctx);
        }
        ctx.f(root, ctx.nodes, ctx.state);
    }

    /**
     * Visit all nodes in a subtree in "post order", meaning leafs get visited first.
     */
    export function doPostOrder<T, S>(tree: ImmutableTree<T>, root: Node<T>, state: S, f: (node: Node<T>, nodes: Nodes<T>, state: S) => boolean | undefined | void) {
        const ctx: VisitorCtx = { nodes: tree.nodes, state, f };
        _doPostOrder(ctx, root);
        return ctx.state;
    }

    function _preOrderFunc(this: VisitorCtx, c: string | undefined) { _doPreOrder(this, this.nodes.get(c!)!); }
    function _doPreOrder<T, S>(ctx: VisitorCtx, root: N) {
        ctx.f(root, ctx.nodes, ctx.state);
        if (root.children.size) {
            root.children.forEach(_preOrderFunc, ctx);
        }
    }

    /**
     * Visit all nodes in a subtree in "pre order", meaning leafs get visited last.
     */
    export function doPreOrder<T, S>(tree: ImmutableTree<T>, root: Node<T>, state: S, f: (node: Node<T>, nodes: Nodes<T>, state: S) => boolean | undefined | void) {
        const ctx: VisitorCtx = { nodes: tree.nodes, state, f };
        _doPreOrder(ctx, root);
        return ctx.state;
    }

    function _subtree(n: N, nodes: Ns, subtree: N[]) { subtree.push(n); }
    /**
     * Get all nodes in a subtree, leafs come first.
     */
    export function subtreePostOrder<T>(tree: ImmutableTree<T>, root: Node<T>) {
        return doPostOrder<T, Node<T>[]>(tree, root, [], _subtree);
    }

    function checkSetRef(oldRef: string, newRef: string) {
        if (oldRef !== newRef) {
            throw new Error(`Cannot setValue of node '${oldRef}' because the new value has a different ref '${newRef}'.`);
        }
    }

    function ensureNotPresent(nodes: Ns, ref: string) {
        if (nodes.has(ref)) {
            throw new Error(`Cannot add node '${ref}' because a different node with this ref already present in the tree.`);
        }
    }

    function ensurePresent(nodes: Ns, ref: string) {
        if (!nodes.has(ref)) {
            throw new Error(`Node '${ref}' is not present in the tree.`);
        }
    }

    function mutateNode(nodes: Ns, mutations: Map<string, N>, ref: string): N {
        ensurePresent(nodes, ref);
        if (mutations.has(ref)) {
            return mutations.get(ref)!;
        }
        const node = nodes.get(ref)!;
        const newNode: N = { ref: node.ref, value: node.value, version: node.version + 1, parent: node.parent, children: node.children.asMutable() };
        mutations.set(ref, newNode);
        nodes.set(ref, newNode);
        return newNode;
    }

    export class Transient<T> implements ImmutableTree<T> {
        nodes = this.tree.nodes.asMutable();
        version: number = this.tree.version + 1;
        private mutations: Map<string, Node<T>> = new Map();

        mutate(ref: string): MutableNode<T> {
            return mutateNode(this.nodes, this.mutations, ref);
        }

        get rootRef() { return this.tree.rootRef; }
        getRef(e: T) {
            return this.tree.getRef(e);
        }

        add(parentRef: string, value: T) {
            const ref = this.getRef(value);
            ensureNotPresent(this.nodes, ref);
            const parent = this.mutate(parentRef);
            const node: Node<T> = { ref, version: 0, value, parent: parent.ref, children: OrderedSet<string>().asMutable() };
            this.mutations.set(ref, node);
            parent.children.add(ref);
            this.nodes.set(ref, node);
            return node;
        }

        setValue(ref: string, value: T): Node<T> {
            checkSetRef(ref, this.getRef(value));
            const node = this.mutate(ref);
            node.value = value;
            return node;
        }

        remove<T>(ref: string): Node<T>[] {
            const { nodes, mutations, mutate } = this;
            const node = nodes.get(ref);
            if (!node) return [];
            const parent = nodes.get(node.parent)!;
            const children = mutate(parent.ref).children;
            const st = subtreePostOrder(this, node);
            if (parent.ref === node.ref) {
                nodes.clear();
                mutations.clear();
                return st;
            }
            children.delete(ref);
            for (const n of st) {
                nodes.delete(n.value.ref);
                mutations.delete(n.value.ref);
            }
            return st;
        }

        removeChildren(ref: string): Node<T>[] {
            const { nodes, mutations, mutate } = this;
            let node = nodes.get(ref);
            if (!node || !node.children.size) return [];
            node = mutate(ref);
            const st = subtreePostOrder(this, node);
            node.children.clear();
            for (const n of st) {
                if (n === node) continue;
                nodes.delete(n.value.ref);
                mutations.delete(n.value.ref);
            }
            return st;
        }

        asImmutable() {
            if (this.mutations.size === 0) return this.tree;

            this.mutations.forEach(m => (m as MutableNode<T>).children = m.children.asImmutable());
            return new Impl<T>(this.tree.rootRef, this.nodes.asImmutable(), this.tree.getRef, this.version);
        }

        constructor(private tree: ImmutableTree<T>) {

        }
    }
}