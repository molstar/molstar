/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Task from '../task'
import RuntimeContext from './runtime-context'

const voidPromise = Promise.resolve(void 0);

class SynchronousRuntimeContext implements RuntimeContext {
    id: number = 0;
    requiresUpdate: boolean = false;
    update(progress: Partial<RuntimeContext.ProgressUpdate>): Promise<void> { return voidPromise; }
    runChild<T>(progress: Partial<RuntimeContext.ProgressUpdate>, task: Task<T>): Promise<T> { return ExecuteSynchronous(task); }
}

const SyncRuntimeInstance = new SynchronousRuntimeContext();

function ExecuteSynchronous<T>(task: Task<T>) {
    return task.__f(SyncRuntimeInstance);
}

export default ExecuteSynchronous