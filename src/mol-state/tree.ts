/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transform } from './transform';
import { ImmutableTree } from './util/immutable-tree';
import { Transformer } from './transformer';
import { StateObject } from './object';

interface StateTree extends ImmutableTree<Transform> { }

namespace StateTree {
    export interface Transient extends ImmutableTree.Transient<Transform> { }
    export interface Serialized extends ImmutableTree.Serialized { }

    function _getRef(t: Transform) { return t.ref; }

    export function create() {
        return ImmutableTree.create<Transform>(Transform.createRoot('<:root:>'), _getRef);
    }

    export function updateParams<T extends Transformer = Transformer>(tree: StateTree, ref: Transform.Ref, params: Transformer.Params<T>): StateTree {
        const t = tree.nodes.get(ref)!.value;
        const newTransform = Transform.updateParams(t, params);
        const newTree = ImmutableTree.asTransient(tree);
        newTree.setValue(ref, newTransform);
        return newTree.asImmutable();
    }

    export function toJSON(tree: StateTree) {
        return ImmutableTree.toJSON(tree, Transform.toJSON) as Serialized;
    }

    export function fromJSON(data: Serialized): StateTree {
        return ImmutableTree.fromJSON(data, _getRef, Transform.fromJSON);
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
            to<A extends StateObject>(ref: Transform.Ref) { return new To<A>(this.state, ref); }
            toRoot<A extends StateObject>() { return new To<A>(this.state, this.state.tree.rootRef as any); }
            delete(ref: Transform.Ref) { this.state.tree.remove(ref); return this; }
            getTree(): StateTree { return this.state.tree.asImmutable(); }
            constructor(tree: StateTree) { this.state = { tree: ImmutableTree.asTransient(tree) } }
        }

        export class To<A extends StateObject> implements Builder {
            apply<T extends Transformer<A, any, any>>(tr: T, params?: Transformer.Params<T>, props?: Partial<Transform.Options>): To<Transformer.To<T>> {
                const t = tr.apply(params, props);
                this.state.tree.add(this.ref, t);
                return new To(this.state, t.ref);
            }

            getTree(): StateTree { return this.state.tree.asImmutable(); }

            constructor(private state: State, private ref: Transform.Ref) {
                if (!this.state.tree.nodes.has(ref)) {
                    throw new Error(`Could not find node '${ref}'.`);
                }
            }
        }
    }
}

export { StateTree }