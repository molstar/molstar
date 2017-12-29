/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../task'
import { RuntimeContext } from './runtime-context'

class SynchronousRuntimeContext implements RuntimeContext {
    shouldUpdate: boolean = false;
    update(progress: string | Partial<RuntimeContext.ProgressUpdate>, dontNotify?: boolean): Promise<void> | void { }
    runChild<T>(task: Task<T>, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T> { return ExecuteSynchronous(task); }
}

const SyncRuntimeInstance = new SynchronousRuntimeContext();

function ExecuteSynchronous<T>(task: Task<T>) {
    return task.__f(SyncRuntimeInstance);
}

export { ExecuteSynchronous }