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
    readonly rootRef: ImmutableTree.Ref,
    readonly version: number,
    readonly nodes: ImmutableTree.Nodes<T>,
    getRef(e: T): ImmutableTree.Ref,
    getValue(ref: ImmutableTree.Ref): T | undefined
}

export namespace ImmutableTree {
    export type Ref = string
    export interface MutableNode<T> { ref: ImmutableTree.Ref, value: T, version: number, parent: ImmutableTree.Ref, children: OrderedSet<ImmutableTree.Ref> }
    export interface Node<T> extends Readonly<MutableNode<T>> { }
    export interface Nodes<T> extends ImmutableMap<ImmutableTree.Ref, Node<T>> { }

    class Impl<T> implements ImmutableTree<T> {
        readonly rootRef: ImmutableTree.Ref;
        readonly version: number;
        readonly nodes: ImmutableTree.Nodes<T>;
        readonly getRef: (e: T) => ImmutableTree.Ref;

        getValue(ref: Ref) {
            const n = this.nodes.get(ref);
            return n ? n.value : void 0;
        }

        constructor(rootRef: ImmutableTree.Ref, nodes: ImmutableTree.Nodes<T>, getRef: (e: T) => ImmutableTree.Ref, version: number) {
            this.rootRef = rootRef;
            this.nodes = nodes;
            this.getRef = getRef;
            this.version = version;
        }
    }

    /**
     * Create an instance of an immutable tree.
     */
    export function create<T>(root: T, getRef: (t: T) => ImmutableTree.Ref): ImmutableTree<T> {
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

    function _postOrderFunc(this: VisitorCtx, c: ImmutableTree.Ref | undefined) { _doPostOrder(this, this.nodes.get(c!)!); }
    function _doPostOrder(ctx: VisitorCtx, root: N) {
        if (root.children.size) {
            root.children.forEach(_postOrderFunc, ctx);
        }
        ctx.f(root, ctx.nodes, ctx.state);
    }

    /**
     * Visit all nodes in a subtree in "post order", meaning leafs get visited first.
     */
    export function doPostOrder<T, S>(tree: ImmutableTree<T>, root: Node<T>, state: S, f: (node: Node<T>, nodes: Nodes<T>, state: S) => boolean | undefined | void): S {
        const ctx: VisitorCtx = { nodes: tree.nodes, state, f };
        _doPostOrder(ctx, root);
        return ctx.state;
    }

    function _preOrderFunc(this: VisitorCtx, c: ImmutableTree.Ref | undefined) { _doPreOrder(this, this.nodes.get(c!)!); }
    function _doPreOrder(ctx: VisitorCtx, root: N) {
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
    export function doPreOrder<T, S>(tree: ImmutableTree<T>, root: Node<T>, state: S, f: (node: Node<T>, nodes: Nodes<T>, state: S) => boolean | undefined | void): S {
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


    function _visitChildToJson(this: Ref[], ref: Ref) { this.push(ref); }
    interface ToJsonCtx { nodes: Ref[], parent: any, children: any, values: any, valueToJSON: (v: any) => any }
    function _visitNodeToJson(this: ToJsonCtx, node: Node<any>) {
        this.nodes.push(node.ref);
        const children: Ref[] = [];
        node.children.forEach(_visitChildToJson as any, children);
        this.parent[node.ref] = node.parent;
        this.children[node.ref] = children;
        this.values[node.ref] = this.valueToJSON(node.value);
    }

    export interface Serialized {
        root: Ref,
        nodes: Ref[],
        parent: { [key: string]: string },
        children: { [key: string]: any },
        values: { [key: string]: any }
    }

    export function toJSON<T>(tree: ImmutableTree<T>, valueToJSON: (v: T) => any): Serialized {
        const ctx: ToJsonCtx = { nodes: [], parent: { }, children: {}, values: {}, valueToJSON };
        tree.nodes.forEach(_visitNodeToJson as any, ctx);
        return {
            root: tree.rootRef,
            nodes: ctx.nodes,
            parent: ctx.parent,
            children: ctx.children,
            values: ctx.values
        };
    }

    export function fromJSON<T>(data: Serialized, getRef: (v: T) => Ref, valueFromJSON: (v: any) => T): ImmutableTree<T> {
        const nodes = ImmutableMap<ImmutableTree.Ref, Node<T>>().asMutable();
        for (const ref of data.nodes) {
            nodes.set(ref, {
                ref,
                value: valueFromJSON(data.values[ref]),
                version: 0,
                parent: data.parent[ref],
                children: OrderedSet(data.children[ref])
            });
        }
        return new Impl(data.root, nodes.asImmutable(), getRef, 0);
    }

    function checkSetRef(oldRef: ImmutableTree.Ref, newRef: ImmutableTree.Ref) {
        if (oldRef !== newRef) {
            throw new Error(`Cannot setValue of node '${oldRef}' because the new value has a different ref '${newRef}'.`);
        }
    }

    function ensureNotPresent(nodes: Ns, ref: ImmutableTree.Ref) {
        if (nodes.has(ref)) {
            throw new Error(`Cannot add node '${ref}' because a different node with this ref already present in the tree.`);
        }
    }

    function ensurePresent(nodes: Ns, ref: ImmutableTree.Ref) {
        if (!nodes.has(ref)) {
            throw new Error(`Node '${ref}' is not present in the tree.`);
        }
    }

    function mutateNode(nodes: Ns, mutations: Map<ImmutableTree.Ref, N>, ref: ImmutableTree.Ref): N {
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
        private mutations: Map<ImmutableTree.Ref, Node<T>> = new Map();

        mutate(ref: ImmutableTree.Ref): MutableNode<T> {
            return mutateNode(this.nodes, this.mutations, ref);
        }

        get rootRef() { return this.tree.rootRef; }
        getRef(e: T) {
            return this.tree.getRef(e);
        }

        getValue(ref: Ref) {
            const n = this.nodes.get(ref);
            return n ? n.value : void 0;
        }

        add(parentRef: ImmutableTree.Ref, value: T) {
            const ref = this.getRef(value);
            ensureNotPresent(this.nodes, ref);
            const parent = this.mutate(parentRef);
            const node: Node<T> = { ref, version: 0, value, parent: parent.ref, children: OrderedSet<string>().asMutable() };
            this.mutations.set(ref, node);
            parent.children.add(ref);
            this.nodes.set(ref, node);
            return node;
        }

        setValue(ref: ImmutableTree.Ref, value: T): Node<T> {
            checkSetRef(ref, this.getRef(value));
            const node = this.mutate(ref);
            node.value = value;
            return node;
        }

        remove(ref: ImmutableTree.Ref): Node<T>[] {
            if (ref === this.rootRef) {
                return this.removeChildren(ref);
            }

            const { nodes, mutations } = this;
            const node = nodes.get(ref);
            if (!node) return [];
            const parent = nodes.get(node.parent)!;
            this.mutate(parent.ref).children.delete(ref);

            const st = subtreePostOrder(this, node);
            for (const n of st) {
                nodes.delete(n.ref);
                mutations.delete(n.ref);
            }

            return st;
        }

        removeChildren(ref: ImmutableTree.Ref): Node<T>[] {
            const { nodes, mutations } = this;
            let node = nodes.get(ref);
            if (!node || !node.children.size) return [];
            node = this.mutate(ref);
            const st = subtreePostOrder(this, node);
            // remove the last node which is the parent
            st.pop();
            node.children.clear();
            for (const n of st) {
                nodes.delete(n.ref);
                mutations.delete(n.ref);
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