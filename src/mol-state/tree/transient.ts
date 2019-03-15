/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';
import { StateTransform } from '../transform';
import { StateTree } from './immutable';
import { StateObjectCell } from 'mol-state/object';
import { shallowEqual } from 'mol-util/object';

export { TransientTree }

class TransientTree implements StateTree {
    transforms = this.tree.transforms as ImmutableMap<StateTransform.Ref, StateTransform>;
    children = this.tree.children as ImmutableMap<StateTransform.Ref, OrderedSet<StateTransform.Ref>>;
    cellStates = this.tree.cellStates as ImmutableMap<StateTransform.Ref, StateObjectCell.State>;

    private changedNodes = false;
    private changedChildren = false;
    private changedStates = false;

    private _childMutations: Map<StateTransform.Ref, OrderedSet<StateTransform.Ref>> | undefined = void 0;

    private get childMutations() {
        if (this._childMutations) return this._childMutations;
        this._childMutations = new Map();
        return this._childMutations;
    }

    private changeStates() {
        if (this.changedStates) return;
        this.changedStates = true;
        this.cellStates = this.cellStates.asMutable();
    }

    private changeNodes() {
        if (this.changedNodes) return;
        this.changedNodes = true;
        this.transforms = this.transforms.asMutable();
    }

    private changeChildren() {
        if (this.changedChildren) return;
        this.changedChildren = true;
        this.children = this.children.asMutable();
    }

    get root() { return this.transforms.get(StateTransform.RootRef)! }

    cellStatesSnapshot() {
        return this.cellStates.asImmutable();
    }

    asTransient() {
        return this.asImmutable().asTransient();
    }

    private addChild(parent: StateTransform.Ref, child: StateTransform.Ref) {
        this.changeChildren();

        if (this.childMutations.has(parent)) {
            this.childMutations.get(parent)!.add(child);
        } else {
            const set = (this.children.get(parent) as OrderedSet<StateTransform.Ref>).asMutable();
            set.add(child);
            this.children.set(parent, set);
            this.childMutations.set(parent, set);
        }
    }

    private removeChild(parent: StateTransform.Ref, child: StateTransform.Ref) {
        this.changeChildren();

        if (this.childMutations.has(parent)) {
            this.childMutations.get(parent)!.remove(child);
        } else {
            const set = (this.children.get(parent) as OrderedSet<StateTransform.Ref>).asMutable();
            set.remove(child);
            this.children.set(parent, set);
            this.childMutations.set(parent, set);
        }
    }

    private clearRoot() {
        const parent = StateTransform.RootRef;
        if (this.children.get(parent).size === 0) return;

        this.changeChildren();

        const set = OrderedSet<StateTransform.Ref>();
        this.children.set(parent, set);
        this.childMutations.set(parent, set);
    }

    changeParent(ref: StateTransform.Ref, newParent: StateTransform.Ref) {
        ensurePresent(this.transforms, ref);

        const old = this.transforms.get(ref);
        this.removeChild(old.parent, ref);
        this.addChild(newParent, ref);
        this.changeNodes();
        this.transforms.set(ref, StateTransform.withParent(old, newParent));
    }

    updateVersion(ref: StateTransform.Ref) {
        ensurePresent(this.transforms, ref);

        const t = this.transforms.get(ref);
        this.changeNodes();
        this.transforms.set(ref, StateTransform.withNewVersion(t));
    }

    add(transform: StateTransform, initialState?: Partial<StateObjectCell.State>) {
        const ref = transform.ref;

        if (this.transforms.has(transform.ref)) {
            const node = this.transforms.get(transform.ref);
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

        this.changeNodes();
        this.transforms.set(ref, transform);

        if (!this.cellStates.has(ref)) {
            this.changeStates();
            if (StateObjectCell.isStateChange(StateObjectCell.DefaultState, initialState)) {
                this.cellStates.set(ref, { ...StateObjectCell.DefaultState, ...initialState });
            } else {
                this.cellStates.set(ref, StateObjectCell.DefaultState);
            }
        }

        return this;
    }

    /** Calls Transform.definition.params.areEqual if available, otherwise uses shallowEqual to check if the params changed */
    setParams(ref: StateTransform.Ref, params: any) {
        ensurePresent(this.transforms, ref);

        const transform = this.transforms.get(ref)!;
        // TODO: should this be here?
        if (shallowEqual(transform.params, params)) {
            return false;
        }

        if (!this.changedNodes) {
            this.changedNodes = true;
            this.transforms = this.transforms.asMutable();
        }

        this.transforms.set(transform.ref, StateTransform.withParams(transform, params));
        return true;
    }

    updateCellState(ref: StateTransform.Ref, state: Partial<StateObjectCell.State>) {
        ensurePresent(this.transforms, ref);

        const old = this.cellStates.get(ref);
        if (!StateObjectCell.isStateChange(old, state)) return false;

        this.changeStates();
        this.cellStates.set(ref, { ...old, ...state });

        return true;
    }

    remove(ref: StateTransform.Ref): StateTransform[] {
        const node = this.transforms.get(ref);
        if (!node) return [];

        const st = StateTree.subtreePostOrder(this, node);
        if (ref === StateTransform.RootRef) {
            st.pop();
            if (st.length === 0) return st;
            this.clearRoot();
        } else {
            if (st.length === 0) return st;
            this.removeChild(node.parent, node.ref);
        }

        this.changeNodes();
        this.changeChildren();
        this.changeStates();

        for (const n of st) {
            this.transforms.delete(n.ref);
            this.children.delete(n.ref);
            this.cellStates.delete(n.ref);
            if (this._childMutations) this._childMutations.delete(n.ref);
        }

        return st;
    }

    asImmutable() {
        if (!this.changedNodes && !this.changedChildren && !this.changedStates && !this._childMutations) return this.tree;
        if (this._childMutations) this._childMutations.forEach(fixChildMutations, this.children);
        return StateTree.create(
            this.changedNodes ? this.transforms.asImmutable() : this.transforms,
            this.changedChildren ? this.children.asImmutable() : this.children,
            this.changedStates ? this.cellStates.asImmutable() : this.cellStates);
    }

    constructor(private tree: StateTree) {

    }
}

function fixChildMutations(this: ImmutableMap<StateTransform.Ref, OrderedSet<StateTransform.Ref>>, m: OrderedSet<StateTransform.Ref>, k: StateTransform.Ref) { this.set(k, m.asImmutable()); }

function alreadyPresent(ref: StateTransform.Ref) {
    throw new Error(`Transform '${ref}' is already present in the tree.`);
}

function parentNotPresent(ref: StateTransform.Ref) {
    throw new Error(`Parent '${ref}' must be present in the tree.`);
}

function ensurePresent(nodes: StateTree.Transforms, ref: StateTransform.Ref) {
    if (!nodes.has(ref)) {
        throw new Error(`Node '${ref}' is not present in the tree.`);
    }
}