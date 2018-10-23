/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ImmutableTree } from '../util/immutable-tree';
import { TransformTree } from './tree';
import { StateObject } from '../object';
import { Transform } from './transform';
import { Transformer } from '../transformer';

export interface StateTreeBuilder {
    getTree(): TransformTree
}

export namespace StateTreeBuilder {
    interface State {
        tree: TransformTree.Transient
    }

    export function create(tree: TransformTree) {
        return new Root(tree);
    }

    export class Root implements StateTreeBuilder {
        private state: State;
        //to<A extends StateObject>(ref: Transform.Ref) { return new To<A>(this.state, ref); }
        to<A extends StateObject>(ref: Transform.Ref) { return new To<A>(this.state, ref); }
        toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.rootRef as any); }
        getTree(): TransformTree { return this.state.tree.asImmutable(); }
        constructor(tree: TransformTree) { this.state = { tree: ImmutableTree.asTransient(tree) } }
    }

    export class To<A extends StateObject> implements StateTreeBuilder {
        apply<T extends Transformer<A, any, any>>(tr: T, params?: Transformer.Params<T>, props?: Partial<Transform.Props>): To<Transformer.To<T>> {
            const t = tr.apply(params, props);
            this.state.tree.add(this.ref, t);
            return new To(this.state, t.ref);
        }

        getTree(): TransformTree { return this.state.tree.asImmutable(); }

        constructor(private state: State, private ref: Transform.Ref) {

        }
    }
}