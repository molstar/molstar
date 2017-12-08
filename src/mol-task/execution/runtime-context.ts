/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Task from '../task'

interface RuntimeContext {
    readonly needsYield: boolean,
    updateProgress(progress: string | Partial<RuntimeContext.ProgressUpdate>): void,
    // Idiomatic usage:
    // if (ctx.needsYield) await ctx.yield({ ... });
    yield(progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<void>,
    // Force the user to pass the progress so that the progress tree can be kept in a "good state".
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

export default RuntimeContext