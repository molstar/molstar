/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transform } from './transform';
import { ImmutableTree } from '../util/immutable-tree';
import { Transformer } from '../transformer';

export interface TransformTree extends ImmutableTree<Transform> { }

export namespace TransformTree {
    export interface Transient extends ImmutableTree.Transient<Transform> { }

    function _getRef(t: Transform) { return t.ref; }

    export function create() {
        return ImmutableTree.create<Transform>(Transform.createRoot('<:root:>'), _getRef);
    }

    export function updateParams<T extends Transformer = Transformer>(tree: TransformTree, ref: Transform.Ref, params: Transformer.Params<T>): TransformTree {
        const t = tree.nodes.get(ref)!.value;
        const newTransform = Transform.updateParams(t, params);
        const newTree = ImmutableTree.asTransient(tree);
        newTree.setValue(ref, newTransform);
        return newTree.asImmutable();
    }

    export function toJSON(tree: TransformTree) {
        return ImmutableTree.toJSON(tree, Transform.toJSON);
    }

    export function fromJSON(data: any): TransformTree {
        return ImmutableTree.fromJSON(data, _getRef, Transform.fromJSON);
    }
}