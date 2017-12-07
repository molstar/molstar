/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Task from '../task'
import RuntimeContext from './runtime-context'

interface Progress {
    taskId: number,
    elapsedMs: { real: number, cpu: number },
    tree: Progress.Node,
    tryAbort?: () => void
}

namespace Progress {
    export interface Node {
        readonly progress: Task.Progress,
        readonly children: ReadonlyArray<Node>
    }

    export interface Observer { (progress: Progress): void }
}

class ObservableExecutor {
    async run<T>(task: Task<T>): Promise<T> {
        const ctx = new ObservableRuntimeContext();
        if (!task.__onAbort) return task.__f(ctx);

        try {
            return await task.__f(ctx);
        } catch (e) {
            if (e === Task.Aborted) task.__onAbort();
            throw e;
        }
    }

    constructor(observer: Progress.Observer, updateRateMs: number) {

    }
}

class ObservableRuntimeContext implements RuntimeContext {
    id: number = 0;
    requiresUpdate: boolean = false;
    update(progress: Partial<RuntimeContext.ProgressUpdate>): Promise<void> { 
        return 0 as any;
    }
    runChild<T>(progress: Partial<RuntimeContext.ProgressUpdate>, task: Task<T>): Promise<T> {
        return 0 as any;
    }
}

function ExecuteObservable<T>(task: Task<T>, observer: Progress.Observer, updateRateMs = 250) {
    return new ObservableExecutor(observer, updateRateMs).run(task);
}

namespace ExecuteObservable {
    export let PRINT_ERRORS_TO_CONSOLE = false;
}

export { ExecuteObservable, Progress }