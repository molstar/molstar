/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'immutable';
import { Transform } from '../transform';
import { StateTree } from './immutable';
import { StateTreeBuilder } from './builder';

export { TransientTree }

class TransientTree implements StateTree {
    nodes = this.tree.nodes;
    children = this.tree.children;

    private changedNodes = false;
    private changedChildren = false;
    private _mutations: Map<Transform.Ref, OrderedSet<Transform.Ref>> = void 0 as any;

    private get mutations() {
        if (this._mutations) return this._mutations;
        this._mutations = new Map();
        return this._mutations;
    }

    get root() { return this.nodes.get(Transform.RootRef)! }

    build(): StateTreeBuilder.Root {
        return new StateTreeBuilder.Root(this);
    }

    asTransient() {
        return this.asImmutable().asTransient();
    }

    private addChild(parent: Transform.Ref, child: Transform.Ref) {
        if (!this.changedChildren) {
            this.changedChildren = true;
            this.children = this.children.asMutable();
        }

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
        if (!this.changedChildren) {
            this.changedChildren = true;
            this.children = this.children.asMutable();
        }

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
        if (!this.changedChildren) {
            this.changedChildren = true;
            this.children = this.children.asMutable();
        }
        this.children.set(parent, set);
        this.mutations.set(parent, set);
    }

    add(transform: Transform) {
        const ref = transform.ref;

        if (this.nodes.has(transform.ref)) {
            const node = this.nodes.get(transform.ref);
            if (node.parent !== transform.parent) alreadyPresent(transform.ref);
        }

        const children = this.children.get(transform.parent);
        if (!children) parentNotPresent(transform.parent);

        if (!children.has(transform.ref)) {
            this.addChild(transform.parent, transform.ref);
        }

        if (!this.children.has(transform.ref)) {
            if (!this.changedChildren) {
                this.changedChildren = true;
                this.children = this.children.asMutable();
            }
            this.children.set(transform.ref, OrderedSet());
        }

        if (!this.changedNodes) {
            this.changedNodes = true;
            this.nodes = this.nodes.asMutable();
        }

        this.nodes.set(ref, transform);
        return this;
    }

    set(transform: Transform) {
        ensurePresent(this.nodes, transform.ref);

        if (!this.changedNodes) {
            this.changedNodes = true;
            this.nodes = this.nodes.asMutable();
        }

        this.nodes.set(transform.ref, transform);
        return this;
    }

    remove(ref: Transform.Ref): Transform[] {
        const node = this.nodes.get(ref);
        if (!node) return [];

        const st = StateTree.subtreePostOrder(this, node);
        if (ref === Transform.RootRef) {
            st.pop();
            if (st.length === 0) return st;
            this.clearRoot();
        } else {
            if (st.length === 0) return st;
            this.removeChild(node.parent, node.ref);
        }

        if (!this.changedNodes) {
            this.changedNodes = true;
            this.nodes = this.nodes.asMutable();
        }

        for (const n of st) {
            this.nodes.delete(n.ref);
            this.children.delete(n.ref);
            if (this._mutations) this._mutations.delete(n.ref);
        }

        return st;
    }

    asImmutable() {
        if (!this.changedNodes && !this.changedChildren && !this._mutations) return this.tree;
        if (this._mutations) this._mutations.forEach(fixChildMutations, this.children);
        return StateTree.create(
            this.changedNodes ? this.nodes.asImmutable() : this.nodes,
            this.changedChildren ? this.children.asImmutable() : this.children);
    }

    constructor(private tree: StateTree) {

    }
}

function fixChildMutations(this: StateTree.Children, m: OrderedSet<Transform.Ref>, k: Transform.Ref) { this.set(k, m.asImmutable()); }

function alreadyPresent(ref: Transform.Ref) {
    throw new Error(`Transform '${ref}' is already present in the tree.`);
}

function parentNotPresent(ref: Transform.Ref) {
    throw new Error(`Parent '${ref}' must be present in the tree.`);
}

function ensurePresent(nodes: StateTree.Nodes, ref: Transform.Ref) {
    if (!nodes.has(ref)) {
        throw new Error(`Node '${ref}' is not present in the tree.`);
    }
}