/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTree } from './immutable';
import { TransientTree } from './transient';
import { StateObject } from '../object';
import { Transform } from '../transform';
import { Transformer } from '../transformer';

export { StateTreeBuilder }

interface StateTreeBuilder {
    readonly editInfo: StateTreeBuilder.EditInfo,
    getTree(): StateTree
}

namespace StateTreeBuilder {
    export interface EditInfo {
        sourceTree: StateTree,
        count: number,
        lastUpdate?: Transform.Ref
    }

    interface State {
        tree: TransientTree,
        editInfo: EditInfo
    }

    export function is(obj: any): obj is StateTreeBuilder {
        return !!obj && typeof (obj as StateTreeBuilder).getTree === 'function';
    }

    export class Root implements StateTreeBuilder {
        private state: State;
        get editInfo() { return this.state.editInfo; }

        to<A extends StateObject>(ref: Transform.Ref) { return new To<A>(this.state, ref, this); }
        toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.root.ref, this); }
        delete(ref: Transform.Ref) {
            this.editInfo.count++;
            this.state.tree.remove(ref);
            return this;
        }
        getTree(): StateTree { return this.state.tree.asImmutable(); }
        constructor(tree: StateTree) { this.state = { tree: tree.asTransient(), editInfo: { sourceTree: tree, count: 0, lastUpdate: void 0 } } }
    }

    export class To<A extends StateObject> implements StateTreeBuilder {
        get editInfo() { return this.state.editInfo; }

        apply<T extends Transformer<A, any, any>>(tr: T, params?: Transformer.Params<T>, props?: Partial<Transform.Options>): To<Transformer.To<T>> {
            const t = tr.apply(this.ref, params, props);
            this.state.tree.add(t);
            this.editInfo.count++;
            this.editInfo.lastUpdate = t.ref;
            return new To(this.state, t.ref, this.root);
        }

        update<T extends Transformer<A, any, any>>(transformer: T, params: (old: Transformer.Params<T>) => Transformer.Params<T>): Root
        update(params: any): Root
        update<T extends Transformer<A, any, any>>(paramsOrTransformer: T, provider?: (old: Transformer.Params<T>) => Transformer.Params<T>) {
            let params: any;
            if (provider) {
                const old = this.state.tree.nodes.get(this.ref)!;
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

        to<A extends StateObject>(ref: Transform.Ref) { return this.root.to<A>(ref); }
        toRoot<A extends StateObject>() { return this.root.toRoot<A>(); }
        delete(ref: Transform.Ref) { return this.root.delete(ref); }

        getTree(): StateTree { return this.state.tree.asImmutable(); }

        constructor(private state: State, private ref: Transform.Ref, private root: Root) {
            if (!this.state.tree.nodes.has(ref)) {
                throw new Error(`Could not find node '${ref}'.`);
            }
        }
    }
}