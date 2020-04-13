/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Map as ImmutableMap, OrderedSet } from 'immutable';
import { StateTransform } from '../transform';
import { StateTree } from './immutable';
import { shallowEqual } from '../../mol-util/object';
import { arrayEqual } from '../../mol-util/array';

export { TransientTree };

class TransientTree implements StateTree {
    transforms = this.tree.transforms as ImmutableMap<StateTransform.Ref, StateTransform>;
    children = this.tree.children as ImmutableMap<StateTransform.Ref, OrderedSet<StateTransform.Ref>>;
    dependencies = this.tree.dependencies as ImmutableMap<StateTransform.Ref, OrderedSet<StateTransform.Ref>>;

    private changedNodes = false;
    private changedChildren = false;
    private changedDependencies = false;

    private _childMutations: Map<StateTransform.Ref, OrderedSet<StateTransform.Ref>> | undefined = void 0;
    private _dependencyMutations: Map<StateTransform.Ref, OrderedSet<StateTransform.Ref>> | undefined = void 0;
    private _stateUpdates: Set<StateTransform.Ref> | undefined = void 0;

    private get childMutations() {
        if (this._childMutations) return this._childMutations;
        this._childMutations = new Map();
        return this._childMutations;
    }

    private get dependencyMutations() {
        if (this._dependencyMutations) return this._dependencyMutations;
        this._dependencyMutations = new Map();
        return this._dependencyMutations;
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

    private changeDependencies() {
        if (this.changedDependencies) return;
        this.changedDependencies = true;
        this.dependencies = this.dependencies.asMutable();
    }

    get root() { return this.transforms.get(StateTransform.RootRef)!; }

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

    private mutateDependency(parent: StateTransform.Ref, child: StateTransform.Ref, action: 'add' | 'remove') {
        let set = this.dependencyMutations.get(parent);

        if (!set) {
            const src = this.dependencies.get(parent);
            if (!src && action === 'remove') return;

            this.changeDependencies();
            set = src ? src.asMutable() : OrderedSet<string>().asMutable();
            this.dependencyMutations.set(parent, set);
            this.dependencies.set(parent, set);
        }

        if (action === 'add') {
            set.add(child);
        } else {
            set.remove(child);
        }
    }

    changeParent(ref: StateTransform.Ref, newParent: StateTransform.Ref) {
        ensurePresent(this.transforms, ref);

        const old = this.transforms.get(ref);
        this.removeChild(old.parent, ref);
        this.addChild(newParent, ref);
        this.changeNodes();
        this.transforms.set(ref, StateTransform.withParent(old, newParent));
    }

    add(transform: StateTransform) {
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

        if (transform.dependsOn) {
            for (const d of transform.dependsOn) {
                this.mutateDependency(d, ref, 'add');
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

    /** Calls Transform.definition.params.areEqual if available, otherwise uses shallowEqual to check if the params changed */
    setTags(ref: StateTransform.Ref, tags: string | string[] | undefined) {
        ensurePresent(this.transforms, ref);

        const transform = this.transforms.get(ref)!;

        const withTags = StateTransform.withParams(transform, tags);
        // TODO: should this be here?
        if (arrayEqual(transform.tags, withTags.tags)) {
            return false;
        }

        if (!this.changedNodes) {
            this.changedNodes = true;
            this.transforms = this.transforms.asMutable();
        }

        this.transforms.set(transform.ref, withTags);
        return true;
    }

    assignState(ref: StateTransform.Ref, state?: Partial<StateTransform.State>) {
        ensurePresent(this.transforms, ref);

        const old = this.transforms.get(ref);
        if (this._stateUpdates && this._stateUpdates.has(ref)) {
            StateTransform.assignState(old.state, state);
            return old;
        } else {
            if (!this._stateUpdates) this._stateUpdates = new Set();
            this._stateUpdates.add(old.ref);
            this.changeNodes();
            const updated = StateTransform.withState(old, state);
            this.transforms.set(ref, updated);
            return updated;
        }
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


        for (const n of st) {
            this.transforms.delete(n.ref);
            this.children.delete(n.ref);
            if (this._childMutations) this._childMutations.delete(n.ref);
        }

        const depRemoves: StateTransform[] = [];
        for (const n of st) {

            if (n.dependsOn) {
                for (const d of n.dependsOn) {
                    if (!this.transforms.has(d)) continue;
                    this.mutateDependency(d, n.ref, 'remove');
                }
            }

            if (this.dependencies.has(n.ref)) {
                const deps = this.dependencies.get(n.ref).toArray();
                this.changeDependencies();
                this.dependencies.delete(n.ref);
                if (this._dependencyMutations) this._dependencyMutations.delete(n.ref);

                for (const dep of deps) {
                    if (!this.transforms.has(dep)) continue;
                    for (const del of this.remove(dep)) depRemoves[depRemoves.length] = del;
                }
            }
        }

        for (const dep of depRemoves) st[st.length] = dep;

        return st;
    }

    asImmutable() {
        if (!this.changedNodes && !this.changedChildren && !this._childMutations) return this.tree;
        if (this._childMutations) this._childMutations.forEach(fixChildMutations, this.children);
        if (this._dependencyMutations) this._dependencyMutations.forEach(fixDependencyMutations, this.dependencies);
        return StateTree.create(
            this.changedNodes ? this.transforms.asImmutable() : this.transforms,
            this.changedChildren ? this.children.asImmutable() : this.children,
            this.changedDependencies ? this.dependencies.asImmutable() : this.dependencies);
    }

    constructor(private tree: StateTree) {

    }
}

function fixChildMutations(this: ImmutableMap<StateTransform.Ref, OrderedSet<StateTransform.Ref>>, m: OrderedSet<StateTransform.Ref>, k: StateTransform.Ref) {
    this.set(k, m.asImmutable());
}

function fixDependencyMutations(this: ImmutableMap<StateTransform.Ref, OrderedSet<StateTransform.Ref>>, m: OrderedSet<StateTransform.Ref>, k: StateTransform.Ref) {
    if (m.size === 0) this.delete(k);
    else this.set(k, m.asImmutable());
}

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