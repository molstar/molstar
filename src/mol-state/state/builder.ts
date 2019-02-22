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

    interface State {
        tree: TransientTree,
        editInfo: EditInfo
    }

    export function is(obj: any): obj is StateBuilder {
        return !!obj && typeof (obj as StateBuilder).getTree === 'function';
    }

    export function isTo(obj: any): obj is StateBuilder.To<any> {
        return !!obj && typeof (obj as StateBuilder).getTree === 'function' && typeof (obj as StateBuilder.To<any>).ref === 'string';
    }

    export class Root implements StateBuilder {
        private state: State;
        get editInfo() { return this.state.editInfo; }

        to<A extends StateObject>(ref: StateTransform.Ref) { return new To<A>(this.state, ref, this); }
        toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.root.ref, this); }
        delete(ref: StateTransform.Ref) {
            this.editInfo.count++;
            this.state.tree.remove(ref);
            return this;
        }
        getTree(): StateTree { return this.state.tree.asImmutable(); }
        constructor(tree: StateTree) { this.state = { tree: tree.asTransient(), editInfo: { sourceTree: tree, count: 0, lastUpdate: void 0 } } }
    }

    export class To<A extends StateObject> implements StateBuilder {
        get editInfo() { return this.state.editInfo; }

        readonly ref: StateTransform.Ref;

        /**
         * Apply the transformed to the parent node
         * If no params are specified (params <- undefined), default params are lazily resolved.
         */
        apply<T extends StateTransformer<A, any, any>>(tr: T, params?: StateTransformer.Params<T>, options?: Partial<StateTransform.Options>, initialCellState?: Partial<StateObjectCell.State>): To<StateTransformer.To<T>> {
            const t = tr.apply(this.ref, params, options);
            this.state.tree.add(t, initialCellState);
            this.editInfo.count++;
            this.editInfo.lastUpdate = t.ref;
            return new To(this.state, t.ref, this.root);
        }

        update<T extends StateTransformer<any, A, any>>(transformer: T, params: (old: StateTransformer.Params<T>) => StateTransformer.Params<T>): Root
        update(params: any): Root
        update<T extends StateTransformer<any, A, any>>(paramsOrTransformer: T, provider?: (old: StateTransformer.Params<T>) => StateTransformer.Params<T>) {
            let params: any;
            if (provider) {
                const old = this.state.tree.transforms.get(this.ref)!;
                params = provider(old.params as any);
            } else {
                params = paramsOrTransformer;
            }

            if (this.state.tree.setParams(this.ref, params)) {
                this.editInfo.count++;
                this.editInfo.lastUpdate = this.ref;
            }

            return this.root;
        }

        to<A extends StateObject>(ref: StateTransform.Ref) { return this.root.to<A>(ref); }
        toRoot<A extends StateObject>() { return this.root.toRoot<A>(); }
        delete(ref: StateTransform.Ref) { return this.root.delete(ref); }

        getTree(): StateTree { return this.state.tree.asImmutable(); }

        constructor(private state: State, ref: StateTransform.Ref, private root: Root) {
            this.ref = ref;
            if (!this.state.tree.transforms.has(ref)) {
                throw new Error(`Could not find node '${ref}'.`);
            }
        }
    }
}