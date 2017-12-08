/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Task from '../task'

interface Progress {
    root: Progress.Node,
    canAbort: boolean,
    tryAbort: (reason?: string) => void
}

namespace Progress {
    export interface Node {
        readonly progress: Task.Progress,
        readonly children: ReadonlyArray<Node>
    }

    export interface Observer { (progress: Progress): void }
}

export default Progress