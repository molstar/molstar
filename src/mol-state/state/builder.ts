/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree } from '../tree/immutable';
import { TransientTree } from '../tree/transient';
import { StateObject, StateObjectCell } from '../object';
import { StateTransform } from '../transform';
import { StateTransformer } from '../transformer';
import { State } from '../state';

export { StateBuilder }

interface StateBuilder {
    readonly editInfo: StateBuilder.EditInfo,
    getTree(): StateTree
}

namespace StateBuilder {
    export interface EditInfo {
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

    export class Root implements StateBuilder {
        private state: BuildState;
        get editInfo() { return this.state.editInfo; }

        to<A extends StateObject>(ref: StateTransform.Ref): To<A>
        to<C extends StateObjectCell>(cell: C): To<StateObjectCell.Obj<C>, StateTransform.Transformer<StateObjectCell.Transform<C>>>
        to(refOrCell: StateTransform.Ref | StateObjectCell) {
            const ref = typeof refOrCell === 'string' ? refOrCell : refOrCell.transform.ref;
            return new To<StateObject, StateTransformer>(this.state, ref, this);
        }
        toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.root.ref, this); }
        delete(ref: StateTransform.Ref) {
            if (!this.state.tree.transforms.has(ref)) return this;
            this.editInfo.count++;
            this.state.tree.remove(ref);
            this.state.actions.push({ kind: 'delete', ref });
            return this;
        }
        getTree(): StateTree { return buildTree(this.state); }
        constructor(tree: StateTree, state?: State) { this.state = { state, tree: tree.asTransient(), actions: [], editInfo: { sourceTree: tree, count: 0, lastUpdate: void 0 } } }
    }

    export class To<A extends StateObject, T extends StateTransformer = StateTransformer> implements StateBuilder {
        get editInfo() { return this.state.editInfo; }

        readonly ref: StateTransform.Ref;

        /**
         * Apply the transformed to the parent node
         * If no params are specified (params <- undefined), default params are lazily resolved.
         */
        apply<T extends StateTransformer<A, any, any>>(tr: T, params?: StateTransformer.Params<T>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>> {
            const t = tr.apply(this.ref, params, options);
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
        applyOrUpdate<T extends StateTransformer<A, any, any>>(ref: StateTransform.Ref, tr: T, params?: StateTransformer.Params<T>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>> {
            if (this.state.tree.transforms.has(ref)) {
                this.to(ref).update(params);
                return this.to(ref) as To<StateTransformer.To<T>>;
            } else {
                return this.apply(tr, params, { ...options, ref });
            }
        }

        /**
         * A helper to greate a group-like state object and keep the current type.
         */
        group<T extends StateTransformer<A, any, any>>(tr: T, params?: StateTransformer.Params<T>, options?: Partial<StateTransform.Options>): To<A> {
            return this.apply(tr, params, options) as To<A>;
        }

        /**
         * Inserts a new transform that does not change the object type and move the original children to it.
         */
        insert<T extends StateTransformer<A, A, any>>(tr: T, params?: StateTransformer.Params<T>, options?: Partial<StateTransform.Options>): To<StateTransformer.To<T>> {
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

        // /**
        //  * Updates a transform in an instantiated tree, passing the transform's source into the providers
        //  *
        //  * This only works if the transform source is NOT updated by the builder. Use at own discression.
        //  */
        // updateInState<T extends StateTransformer<any, A, any>>(transformer: T, state: State, provider: (old: StateTransformer.Params<T>, a: StateTransformer.From<T>) => StateTransformer.Params<T>): Root {
        //     const old = this.state.tree.transforms.get(this.ref)!;
        //     const cell = state.cells.get(this.ref);
        //     if (!cell || !cell.sourceRef) throw new Error('Source cell is not present in the tree.');
        //     const parent = state.cells.get(cell.sourceRef);
        //     if (!parent || !parent.obj) throw new Error('Parent cell is not present or computed.');

        //     const params = provider(old.params as any, parent.obj as any);

        //     if (this.state.tree.setParams(this.ref, params)) {
        //         this.editInfo.count++;
        //         this.editInfo.lastUpdate = this.ref;
        //     }

        //     return this.root;
        // }

        update<T extends StateTransformer<any, A, any>>(transformer: T, params: (old: StateTransformer.Params<T>) => StateTransformer.Params<T>): Root
        update(params: StateTransformer.Params<T> | ((old: StateTransformer.Params<T>) => StateTransformer.Params<T>)): Root
        update<T extends StateTransformer<any, A, any>>(paramsOrTransformer: T | any, provider?: (old: StateTransformer.Params<T>) => StateTransformer.Params<T>) {
            let params: any;
            if (provider) {
                const old = this.state.tree.transforms.get(this.ref)!;
                params = provider(old.params as any);
            } else {
                params = typeof paramsOrTransformer === 'function'
                    ? paramsOrTransformer(this.state.tree.transforms.get(this.ref)!.params)
                    : paramsOrTransformer;
            }

            if (this.state.tree.setParams(this.ref, params)) {
                this.editInfo.count++;
                this.editInfo.lastUpdate = this.ref;
            }

            this.state.actions.push({ kind: 'update', ref: this.ref, params });

            return this.root;
        }

        to<A extends StateObject>(ref: StateTransform.Ref) { return this.root.to<A>(ref); }
        toRoot<A extends StateObject>() { return this.root.toRoot<A>(); }
        delete(ref: StateTransform.Ref) { return this.root.delete(ref); }

        getTree(): StateTree { return buildTree(this.state); }

        constructor(private state: BuildState, ref: StateTransform.Ref, private root: Root) {
            this.ref = ref;
            if (!this.state.tree.transforms.has(ref)) {
                throw new Error(`Could not find node '${ref}'.`);
            }
        }
    }
}