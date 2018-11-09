/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transform } from './transform';
import { ImmutableTree, TransientTree } from './tree/base';
import { Transformer } from './transformer';
import { StateObject } from './object';

export { StateTree, ImmutableTree, TransientTree }

interface StateTree extends ImmutableTree { }

namespace StateTree {
    export interface Transient extends TransientTree { }
    export interface Serialized extends ImmutableTree.Serialized { }

    export function create() {
        return ImmutableTree.create(Transform.createRoot(Transform.RootRef));
    }

    export function updateParams<T extends Transformer = Transformer>(tree: StateTree, ref: Transform.Ref, params: Transformer.Params<T>): StateTree {
        const t = tree.nodes.get(ref)!;
        const newTransform = Transform.updateParams(t, params);
        return tree.asTransient().set(newTransform).asImmutable();
    }

    export function toJSON(tree: StateTree) {
        return ImmutableTree.toJSON(tree) as Serialized;
    }

    export function fromJSON(data: Serialized): StateTree {
        return ImmutableTree.fromJSON(data);
    }

    export interface Builder {
        getTree(): StateTree
    }

    export function build(tree: StateTree) {
        return new Builder.Root(tree);
    }

    export namespace Builder {
        interface State {
            tree: StateTree.Transient
        }

        export class Root implements Builder {
            private state: State;
            to<A extends StateObject>(ref: Transform.Ref) { return new To<A>(this.state, ref, this); }
            toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.root.ref, this); }
            delete(ref: Transform.Ref) {
                this.state.tree.remove(ref);
                return this;
            }
            getTree(): StateTree { return this.state.tree.asImmutable(); }
            constructor(tree: StateTree) { this.state = { tree: tree.asTransient() } }
        }

        export class To<A extends StateObject> implements Builder {
            apply<T extends Transformer<A, any, any>>(tr: T, params?: Transformer.Params<T>, props?: Partial<Transform.Options>): To<Transformer.To<T>> {
                const t = tr.apply(this.ref, params, props);
                this.state.tree.add(t);
                return new To(this.state, t.ref, this.root);
            }

            and() { return this.root; }

            getTree(): StateTree { return this.state.tree.asImmutable(); }

            constructor(private state: State, private ref: Transform.Ref, private root: Root) {
                if (!this.state.tree.nodes.has(ref)) {
                    throw new Error(`Could not find node '${ref}'.`);
                }
            }
        }
    }
}