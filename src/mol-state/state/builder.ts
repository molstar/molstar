/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree } from '../tree/immutable';
import { TransientTree } from '../tree/transient';
import { StateObject, StateObjectCell, StateObjectSelector, StateObjectRef } from '../object';
import { StateTransform } from '../transform';
import { StateTransformer } from '../transformer';
import { State } from '../state';
import { produce } from 'immer';

export { StateBuilder };

interface StateBuilder {
    readonly editInfo: StateBuilder.EditInfo,
    getTree(): StateTree
}

namespace StateBuilder {
    export interface EditInfo {
        applied: boolean,
        sourceTree: StateTree,
        count: number,
        lastUpdate?: StateTransform.Ref
    }

    interface BuildState {
        state: State | undefined,
        tree: TransientTree,
        editInfo: EditInfo,
        actions: Action[]
    }

    type Action =
        | { kind: 'add', transform: StateTransform }
        | { kind: 'update', ref: string, params: any }
        | { kind: 'delete', ref: string }
        | { kind: 'insert', ref: string, transform: StateTransform }

    function buildTree(state: BuildState) {
        if (!state.state || state.state.tree === state.editInfo.sourceTree) {
            return state.tree.asImmutable();
        }

        // The tree has changed in the meantime, we need to reapply the changes!
        const tree = state.state.tree.asTransient();
        for (const a of state.actions) {
            switch (a.kind) {
                case 'add': tree.add(a.transform); break;
                case 'update': tree.setParams(a.ref, a.params); break;
                case 'delete': tree.remove(a.ref); break;
                case 'insert': {
                    const children = tree.children.get(a.ref).toArray();
                    tree.add(a.transform);
                    for (const c of children) {
                        tree.changeParent(c, a.transform.ref);
                    }
                    break;
                }
            }
        }
        state.editInfo.sourceTree = state.tree;
        return tree.asImmutable();
    }

    export function is(obj: any): obj is StateBuilder {
        return !!obj && typeof (obj as StateBuilder).getTree === 'function';
    }

    export function isTo(obj: any): obj is StateBuilder.To<any> {
        return !!obj && typeof (obj as StateBuilder).getTree === 'function' && typeof (obj as StateBuilder.To<any>).ref === 'string';
    }

    // type ToFromCell<C extends StateObjectCell> = C extends StateObjectCell<infer A, StateTransform<infer T extends StateTransformer>> ? To<A, any>: never

    export class Root implements StateBuilder {
        private state: BuildState;
        get editInfo() { return this.state.editInfo; }
        get currentTree() { return this.state.tree; }

        to<A extends StateObject, T extends StateTransformer>(ref: StateTransform.Ref): To<A, T>
        to<A extends StateObject, T extends StateTransformer>(ref: StateObjectRef<A>): To<A, T>
        to<C extends StateObjectCell>(cell: C): To<StateObjectCell.Obj<C>, StateObjectCell.Transformer<C>>
        to<S extends StateObjectSelector>(selector: S): To<StateObjectSelector.Obj<S>, StateObjectSelector.Transformer<S>>
        to(refOrCellOrSelector: StateTransform.Ref | StateObjectCell | StateObjectSelector) {
            const ref = typeof refOrCellOrSelector === 'string'
                ? refOrCellOrSelector
                : StateObjectCell.is(refOrCellOrSelector)
                    ? refOrCellOrSelector.transform.ref
                    : refOrCellOrSelector.ref;
            return new To<StateObject, StateTransformer>(this.state, ref, this);
        }
        toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.root.ref, this); }
        delete(obj: StateObjectRef) {
            const ref = StateObjectRef.resolveRef(obj);
            if (!ref || !this.state.tree.transforms.has(ref)) return this;
            this.editInfo.count++;
            this.state.tree.remove(ref);
            this.state.actions.push({ kind: 'delete', ref });
            return this;
        }
        getTree(): StateTree { return buildTree(this.state); }

        commit(options?: Partial<State.UpdateOptions>) {
            if (!this.state.state) throw new Error('Cannot commit template tree');
            return this.state.state.runTask(this.state.state.updateTree(this, options));
        }

        constructor(tree: StateTree, state?: State) { this.state = { state, tree: tree.asTransient(), actions: [], editInfo: { applied: false, sourceTree: tree, count: 0, lastUpdate: void 0 } }; }
    }

    export class To<A extends StateObject, T extends StateTransformer = StateTransformer> implements StateBuilder {
        get editInfo() { return this.state.editInfo; }
        get selector() { return new StateObjectSelector<A, T>(this.ref, this.state.state); }

        readonly ref: StateTransform.Ref;

        private getApplyRoot(): StateTransform.Ref {
            return StateTree.getDecoratorRoot(this.state.tree, this.ref);
        }

        /**
         * Apply the transformed to the parent node
         * If no params are specified (params <- undefined), default params are lazily resolved.
         */
        apply<T extends StateTransformer<A, any, any>>(tr: T, params?: Partial<StateTransformer.Params<T>>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>, T> {
            if (tr.definition.isDecorator) {
                return this.insert(tr, params, options);
            }

            const applyRoot = this.getApplyRoot();
            const t = tr.apply(applyRoot, params, options);

            this.state.tree.add(t);
            this.editInfo.count++;
            this.editInfo.lastUpdate = t.ref;

            this.state.actions.push({ kind: 'add', transform: t });
            return new To(this.state, t.ref, this.root);
        }

        /**
         * If the ref is present, the transform is applied.
         * Otherwise a transform with the specifed ref is created.
         */
        applyOrUpdate<T extends StateTransformer<A, any, any>>(ref: StateTransform.Ref, tr: T, params?: Partial<StateTransformer.Params<T>>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>, T> {
            if (this.state.tree.transforms.has(ref)) {
                const to = this.to<StateTransformer.To<T>, T>(ref);
                if (params) to.update(params);
                return to;
            } else {
                return this.apply(tr, params, { ...options, ref });
            }
        }

        /**
         * Apply the transformed to the parent node
         * If no params are specified (params <- undefined), default params are lazily resolved.
         * The transformer cannot be a decorator to be able to use this.
         */
        applyOrUpdateTagged<T extends StateTransformer<A, any, any>>(tags: string | string[], tr: T, params?: Partial<StateTransformer.Params<T>>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>, T> {
            if (tr.definition.isDecorator) {
                throw new Error(`Can't use applyOrUpdateTagged on decorator transformers.`);
            }

            const applyRoot = this.getApplyRoot();

            const children = this.state.tree.children.get(applyRoot).values();
            while (true) {
                const child = children.next();
                if (child.done) break;
                const tr = this.state.tree.transforms.get(child.value);
                if (tr && StateTransform.hasTags(tr, tags)) {
                    const to = this.to<StateTransformer.To<T>, T>(child.value);
                    to.updateTagged(params, tagsUnion(tr.tags, tags, options && options.tags));
                    return to;
                }
            }

            const t = tr.apply(applyRoot, params, { ...options, tags: tagsUnion(tags, options && options.tags) });
            this.state.tree.add(t);
            this.editInfo.count++;
            this.editInfo.lastUpdate = t.ref;

            this.state.actions.push({ kind: 'add', transform: t });

            return new To(this.state, t.ref, this.root);
        }

        /**
         * A helper to greate a group-like state object and keep the current type.
         */
        group<T extends StateTransformer<A, any, any>>(tr: T, params?: StateTransformer.Params<T>, options?: Partial<StateTransform.Options>): To<A, T> {
            return this.apply(tr, params, options);
        }

        /**
         * Inserts a new transform that does not change the object type and move the original children to it.
         */
        insert<T extends StateTransformer<A, A, any>>(tr: T, params?: Partial<StateTransformer.Params<T>>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>, T> {
            // cache the children
            const children = this.state.tree.children.get(this.ref).toArray();

            // add the new node
            const t = tr.apply(this.ref, params, options);
            this.state.tree.add(t);

            // move the original children to the new node
            for (const c of children) {
                this.state.tree.changeParent(c, t.ref);
            }

            this.editInfo.count++;
            this.editInfo.lastUpdate = t.ref;

            this.state.actions.push({ kind: 'insert', ref: this.ref, transform: t });

            return new To(this.state, t.ref, this.root);
        }

        private updateTagged(params: any, tags: string | string[] | undefined) {
            if (this.state.tree.setParams(this.ref, params) || this.state.tree.setTags(this.ref, tags)) {
                this.editInfo.count++;
                this.editInfo.lastUpdate = this.ref;
                this.state.actions.push({ kind: 'update', ref: this.ref, params });
            }
        }

        update<T extends StateTransformer<any, A, any>>(transformer: T, params: (old: StateTransformer.Params<T>) => Partial<StateTransformer.Params<T>> | void): Root
        update(params: Partial<StateTransformer.Params<T>> | ((old: StateTransformer.Params<T>) => Partial<StateTransformer.Params<T>> | void)): Root
        update<T extends StateTransformer<any, A, any>>(paramsOrTransformer: T | any, provider?: (old: StateTransformer.Params<T>) => StateTransformer.Params<T>) {
            let params: any;
            if (provider) {
                const old = this.state.tree.transforms.get(this.ref)!;
                params = produce(old.params, provider);
            } else {
                params = typeof paramsOrTransformer === 'function'
                    ? produce(this.state.tree.transforms.get(this.ref)!.params, paramsOrTransformer)
                    : paramsOrTransformer;
            }

            if (this.state.tree.setParams(this.ref, params)) {
                this.editInfo.count++;
                this.editInfo.lastUpdate = this.ref;
                this.state.actions.push({ kind: 'update', ref: this.ref, params });
            }
            return this.root;
        }

        to<A extends StateObject, T extends StateTransformer>(ref: StateTransform.Ref): To<A, T>
        to<C extends StateObjectCell>(cell: C): To<StateObjectCell.Obj<C>, StateObjectCell.Transformer<C>>
        to<S extends StateObjectSelector>(selector: S): To<StateObjectSelector.Obj<S>, StateObjectSelector.Transformer<S>>
        to(ref: StateTransform.Ref | StateObjectCell | StateObjectSelector) { return  this.root.to(ref as any); }
        toRoot<A extends StateObject>() { return this.root.toRoot<A>(); }
        delete(ref: StateObjectRef) { return this.root.delete(ref); }

        getTree(): StateTree { return buildTree(this.state); }

        /** Returns selector to this node. */
        commit(options?: Partial<State.UpdateOptions>): Promise<StateObjectSelector<A>> {
            if (!this.state.state) throw new Error('Cannot commit template tree');
            return this.state.state.runTask(this.state.state.updateTree(this, options));
        }

        constructor(private state: BuildState, ref: StateTransform.Ref, private root: Root) {
            this.ref = ref;
            if (!this.state.tree.transforms.has(ref)) {
                throw new Error(`Could not find node '${ref}'.`);
            }
        }
    }
}

function tagsUnion(...arrays: (string[] | string | undefined)[]): string[] | undefined {
    let set: Set<string> | undefined = void 0;
    const ret = [];
    for (const xs of arrays) {
        if (!xs) continue;
        if (!set) set = new Set();
        if (typeof xs === 'string') {
            if (set.has(xs)) continue;
            set.add(xs);
            ret.push(xs);
        } else {
            for (const x of xs) {
                if (set.has(x)) continue;
                set.add(x);
                ret.push(x);
            }
        }
    }
    return ret;
}