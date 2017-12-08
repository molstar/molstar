/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Task from '../task'
import RuntimeContext from './runtime-context'
import Progress from './progress'
import now from '../util/now'

function defaultProgress(rootTaskId: number, task: Task<any>): Task.Progress {
    return {
        rootTaskId,
        taskId: task.id,
        taskName: task.name,
        message: 'Running...',
        elapsedMs: { real: 0, cpu: 0 },
        canAbort: true,
        isIndeterminate: true,
        current: 0,
        max: 0
    };
}

class ProgressInfo {
    taskId: number;
    elapsedMs: { real: number, cpu: number };
    tree: Progress.Node;
    tryAbort?: (reason?: string) => void;

    snapshot(): Progress {
        return 0 as any;
    }
}

class ObservableExecutor {
    progressInfo: ProgressInfo;

    async run<T>(task: Task<T>): Promise<T> {
        const ctx = new ObservableRuntimeContext(task.id, task, 0);
        if (!task.__onAbort) return task.__f(ctx);

        try {
            return await task.__f(ctx);
        } catch (e) {
            if (Task.isAborted(e)) task.__onAbort();
            throw e;
        }
    }

    constructor(observer: Progress.Observer, updateRateMs: number) {

    }
}

class ObservableRuntimeContext implements RuntimeContext {
    elapsedCpuMs: number;
    lastScheduledTime: number;

    started: number;
    taskId: number;
    taskName: string;
    progress: Task.Progress;
    updateRateMs: number;

    get requiresUpdate(): boolean {
        return now() - this.started > this.updateRateMs;
    }

    update(progress: Partial<RuntimeContext.ProgressUpdate>): Promise<void> {
        return 0 as any;
    }

    runChild<T>(progress: Partial<RuntimeContext.ProgressUpdate>, task: Task<T>): Promise<T> {
        return 0 as any;
    }

    constructor(parentId: number, task: Task<any>, updateRateMs: number) {
        this.started = now();
        this.lastScheduledTime = this.started;
        this.taskId = task.id;
        this.taskName = task.name;
        this.progress = defaultProgress(parentId, task);
        this.updateRateMs = updateRateMs;
    }
}

function ExecuteObservable<T>(task: Task<T>, observer: Progress.Observer, updateRateMs = 250) {
    return new ObservableExecutor(observer, updateRateMs).run(task);
}

namespace ExecuteObservable {
    export let PRINT_ERRORS_TO_CONSOLE = false;
}

export default ExecuteObservable