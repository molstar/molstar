/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { BehaviorSubject } from 'rxjs';
import { UUID } from 'mol-util'
import { AnyEntity } from './entity';
import Viewer from '../viewer';
import { Progress } from 'mol-task';

// TODO
export type StateTree = {}

export class StateContext {
    id = UUID.create()
    change = new BehaviorSubject(0)

    tree: StateTree = {}
    entities: Set<AnyEntity> = new Set()

    viewer: Viewer

    constructor(readonly log: (p: Progress) => void) {

    }
}
