/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { TransformTree } from './tree/tree';
import { Transform } from './tree/transform';
import { Map as ImmutableMap } from 'immutable';
import { StateContext } from './context/context';

export interface State<ObjectProps = unknown> {
    definition: State.Definition<ObjectProps>,
    objects: Map<Transform.InstanceId, StateObject>
}

export namespace State {
    export type ObjectProps<P> = ImmutableMap<Transform.InstanceId, P>

    export interface Definition<P = unknown> {
        tree: TransformTree,
        // things like object visibility
        props: ObjectProps<P>
    }

    export async function update<P>(context: StateContext, old: State<P>, tree: Definition<P>, props?: ObjectProps<P>): Promise<State<P>> {
        throw 'nyi';
    }
}
