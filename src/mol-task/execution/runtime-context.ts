/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from './synchronous';

interface RuntimeContext {
    readonly shouldUpdate: boolean,
    readonly isSynchronous: boolean,

    // Idiomatic usage:
    // if (ctx.needsYield) await ctx.yield({ ... });
    //
    // Alternatively, progress can be updated without notifying (and yielding) using update(progress, true).
    // This is useful for nested tasks.
    update(progress?: string | Partial<RuntimeContext.ProgressUpdate>, dontNotify?: boolean): Promise<void> | void
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

    export const Synchronous = SyncRuntimeContext;
}

export { RuntimeContext };