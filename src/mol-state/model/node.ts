/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface ModelNode<T = any> {
    '@type': T
}

export namespace ModelNode {
    export type TypeOf<T>
        = T extends ModelNode<infer X> ? [X]
        : T extends [ModelNode<infer X>] ? [X]
        : T extends [ModelNode<infer X>, ModelNode<infer Y>] ? [X, Y]
        : unknown[];
}