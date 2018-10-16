/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface TransformTree {
    // TODO
}

export namespace TransformTree {
    export interface Update {
        readonly tree: TransformTree,
        readonly rootId: number,
        readonly params: unknown
    }
}