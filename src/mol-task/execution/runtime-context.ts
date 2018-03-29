/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../task'

interface RuntimeContext {
    readonly shouldUpdate: boolean,
    readonly isSynchronous: boolean,

    // Idiomatic usage:
    // if (ctx.needsYield) await ctx.yield({ ... });
    //
    // Alternatively, progress can be updated without notifying (and yielding) using update(progress, true).
    // This is useful for nested tasks.
    update(progress?: string | Partial<RuntimeContext.ProgressUpdate>, dontNotify?: boolean): Promise<void> | void,

    // Run a child task that adds a new node to the progress tree.
    // Allow to pass the progress so that the progress tree can be kept in a "good state" without having to separately call update.
    runChild<T>(task: Task<T>, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T>
}

namespace RuntimeContext {
    export interface AbortToken { isAborted: boolean }

    export interface ProgressUpdate {
        message: string,
        isIndeterminate: boolean,
        current: number,
        max: number,
        canAbort: boolean
    }
}

export { RuntimeContext }