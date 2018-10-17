/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateObject } from './object';
import { TransformTree } from './tree/tree';
import { Transform } from './tree/transform';

export interface State {
    tree: TransformTree,
    objects: Map<Transform.InstanceId, StateObject>,
    history: TransformTree[]
}