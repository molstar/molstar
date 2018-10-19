/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EventDispatcher } from '../event/event';

export interface TransformContext {
    /** An event dispatcher for executing child tasks. */
    dispatcher: EventDispatcher,

    globalContext: any
    // tree: ModelTree
}