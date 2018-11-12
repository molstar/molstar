/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';
import { Transform } from '../transform';
import { StateTree } from './immutable';
import { StateTreeBuilder } from './builder';
import { StateObjectCell } from 'mol-state/object';
import { shallowEqual } from 'mol-util/object';
import { UUID } from 'mol-util';

export { TransientTree }

class TransientTree implements StateTree {
    nodes = this.tree.nodes as ImmutableMap<Transform.Ref, Transform>;
    children = this.tree.children as ImmutableMap<Transform.Ref, OrderedSet<Transform.Ref>>;

    private changedNodes = false;
    private changedChildren = false;
    private _childMutations: Map<Transform.Ref, OrderedSet<Transform.Ref>> | undefined = void 0;
    private _transformMutations: Map<Transform.Ref, Transform> | undefined = void 0;

    private get childMutations() {
        if (this._childMutations) return this._childMutations;
        this._childMutations = new Map();
        return this._childMutations;
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

        if (this.childMutations.has(parent)) {
            this.childMutations.get(parent)!.add(child);
        } else {
            const set = (this.children.get(parent) as OrderedSet<Transform.Ref>).asMutable();
            set.add(child);
            this.children.set(parent, set);
            this.childMutations.set(parent, set);
        }
    }

    private removeChild(parent: Transform.Ref, child: Transform.Ref) {
        if (!this.changedChildren) {
            this.changedChildren = true;
            this.children = this.children.asMutable();
        }

        if (this.childMutations.has(parent)) {
            this.childMutations.get(parent)!.remove(child);
        } else {
            const set = (this.children.get(parent) as OrderedSet<Transform.Ref>).asMutable();
            set.remove(child);
            this.children.set(parent, set);
            this.childMutations.set(parent, set);
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
        this.childMutations.set(parent, set);
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

    /** Calls Transform.definition.params.areEqual if available, otherwise uses shallowEqual to check if the params changed */
    setParams(ref: Transform.Ref, params: unknown) {
        ensurePresent(this.nodes, ref);

        const transform = this.nodes.get(ref)!;
        const def = transform.transformer.definition;
        if (def.params && def.params.areEqual) {
            if (def.params.areEqual(transform.params, params)) return false;
        } else {
            if (shallowEqual(transform.params, params)) {
                return false;
            }
        }

        if (this._transformMutations && this._transformMutations.has(transform.ref)) {
            const mutated = this._transformMutations.get(transform.ref)!;
            (mutated.params as any) = params;
            (mutated.version as UUID) = UUID.create22();
        } else {
            this.set(Transform.withParams(transform, params));
        }

        return true;
    }

    setCellState(ref: Transform.Ref, state: Partial<StateObjectCell.State>) {
        ensurePresent(this.nodes, ref);

        if (this._transformMutations && this._transformMutations.has(ref)) {
            const transform = this._transformMutations.get(ref)!;
            const old = transform.cellState;
            (transform.cellState as StateObjectCell.State) = { ...old, ...state };
            return transform;
        } else {
            const transform = this.nodes.get(ref);
            const newT = Transform.withCellState(transform, state);
            this.set(newT);
            return newT;
        }
    }

    private set(transform: Transform) {
        ensurePresent(this.nodes, transform.ref);

        if (!this.changedNodes) {
            this.changedNodes = true;
            this.nodes = this.nodes.asMutable();
        }

        if (!this._transformMutations) {
            this._transformMutations = new Map();
        }
        this._transformMutations.set(transform.ref, transform);

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
            if (this._childMutations) this._childMutations.delete(n.ref);
        }

        return st;
    }

    asImmutable() {
        if (!this.changedNodes && !this.changedChildren && !this._childMutations) return this.tree;
        if (this._childMutations) this._childMutations.forEach(fixChildMutations, this.children);
        return StateTree.create(
            this.changedNodes ? this.nodes.asImmutable() : this.nodes,
            this.changedChildren ? this.children.asImmutable() : this.children);
    }

    constructor(private tree: StateTree) {

    }
}

function fixChildMutations(this: ImmutableMap<Transform.Ref, OrderedSet<Transform.Ref>>, m: OrderedSet<Transform.Ref>, k: Transform.Ref) { this.set(k, m.asImmutable()); }

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