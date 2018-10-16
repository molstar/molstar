/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Transform } from './transform';
import { ModelNode } from '../model/node';
import { Transformer } from './transformer';

export interface Transform<A extends ModelNode, B extends ModelNode, P = any> {
    readonly instanceId: number,

    readonly transformer: Transformer<A, B, P>,
    readonly props: Transform.Props,

    readonly transformerId: string,
    readonly params: P,
    readonly ref: string,
    readonly version: number,
}

export namespace Transform {
    export interface Props {

    }
}