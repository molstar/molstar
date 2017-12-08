/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Task from '../task'
import RuntimeContext from './runtime-context'
import Progress from './progress'
import now from '../util/now'
import ImmediateScheduler from '../scheduler/immediate'

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

interface ProgressInfo {
    updateRateMs: number,
    lastUpdated: number,
    observer: Progress.Observer,

    abortToken: { abortRequested: boolean, reason: string },
    taskId: number;
    root: Progress.Node;
    tryAbort: (reason?: string) => void;
}

function ProgressInfo(task: Task<any>, observer: Progress.Observer, updateRateMs: number): ProgressInfo {
    const abortToken: ProgressInfo['abortToken'] = { abortRequested: false, reason: '' };

    return {
        updateRateMs,
        lastUpdated: now(),
        observer,
        abortToken,
        taskId: task.id,
        root: { progress: defaultProgress(task.id, task), children: [] },
        tryAbort: abortFn(abortToken)
    };
}

function abortFn(token: ProgressInfo['abortToken']) {
    return (reason?: string) => {
        token.abortRequested = true;
        token.reason = reason || token.reason;
    };
}

function cloneTree(root: Progress.Node): Progress.Node {
    return { progress: { ...root.progress, elapsedMs: { ...root.progress.elapsedMs } }, children: root.children.map(cloneTree) };
}

function canAbort(root: Progress.Node): boolean {
    return root.progress.canAbort && root.children.every(canAbort);
}

function snapshotProgress(info: ProgressInfo): Progress {
    return { root: cloneTree(info.root), canAbort: canAbort(info.root), tryAbort: info.tryAbort };
}


async function runRoot<T>(task: Task<T>, ctx: ObservableRuntimeContext) {
    ctx.started = now();
    if (!task.__onAbort) return task.__f(ctx);
    try {
        const ret = await task.__f(ctx);
        // if (ctx.info.abortToken.abortRequested) {
        //     task.__onAbort();
        // }
        return ret;
    } catch (e) {
        if (Task.isAborted(e)) task.__onAbort();
        if (ExecuteObservable.PRINT_ERRORS_TO_STD_ERR) console.error(e);
        throw e;
    }
}

async function run<T>(task: Task<T>, ctx: ObservableRuntimeContext) {
    ctx.started = now();
    if (!task.__onAbort) return task.__f(ctx);
    try {
        const ret = await task.__f(ctx);
        // if (ctx.info.abortToken.abortRequested) {
        //     task.__onAbort();
        // }
        return ret;
    } catch (e) {
        if (Task.isAborted(e)) task.__onAbort();
        if (ExecuteObservable.PRINT_ERRORS_TO_STD_ERR) console.error(e);
        throw e;
    }
}

class ObservableRuntimeContext implements RuntimeContext {
    isExecuting = true;
    elapsedCpuMs: number;
    lastScheduledTime: number;

    started: number = 0;
    node: Progress.Node;
    info: ProgressInfo;

    private checkAborted() {
        if (this.info.abortToken.abortRequested) {
            throw Task.Aborted(this.info.abortToken.reason);
        }
    }

    get needsYield(): boolean {
        this.checkAborted();
        return now() - this.info.lastUpdated > this.info.updateRateMs;
    }

    private setProgress(update?: string | Partial<RuntimeContext.ProgressUpdate>) {
        this.checkAborted();

        if (!update) return;

        const progress = this.node.progress;
        if (typeof update === 'string') {
            progress.message = update;
        } else {
            if (typeof update.canAbort !== 'undefined') progress.canAbort = update.canAbort;
            if (typeof update.current !== 'undefined') progress.current = update.current;
            if (typeof update.isIndeterminate !== 'undefined') progress.isIndeterminate = update.isIndeterminate;
            if (typeof update.max !== 'undefined') progress.max = update.max;
            if (typeof update.message !== 'undefined') progress.message = update.message;
        }
    }

    private resume = () => {
        this.isExecuting = true;
        this.lastScheduledTime = now();
    }

    updateProgress(progress?: string | Partial<RuntimeContext.ProgressUpdate>) {
        this.setProgress(progress);
    }

    yield(progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<void> {
        this.isExecuting = false;
        this.setProgress(progress);
        this.info.lastUpdated = now();
        const snapshot = snapshotProgress(this.info);
        this.info.observer(snapshot);
        return ImmediateScheduler.last(this.resume);
    }

    async runChild<T>(task: Task<T>, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T> {
        this.setProgress(progress);
        const node: Progress.Node = { progress: defaultProgress(this.info.taskId, task), children: [] };
        const children = this.node.children as Progress.Node[];
        children.push(node);
        const ctx = new ObservableRuntimeContext(task, this.info, node);
        try {
            return await run(task, ctx);
        } finally {
            // remove the progress node after the computation has finished.
            const idx = children.indexOf(node);
            if (idx >= 0) {
                children[idx] = children[children.length - 1];
                children.pop();
            }
        }
    }

    constructor(task: Task<any>, info: ProgressInfo, node: Progress.Node) {
        this.lastScheduledTime = this.started;
        this.node = node;
        this.info = info;
    }
}

function ExecuteObservable<T>(task: Task<T>, observer: Progress.Observer, updateRateMs = 250) {
    const info = ProgressInfo(task, observer, updateRateMs);
    const ctx = new ObservableRuntimeContext(task, info, info.root);
    return runRoot(task, ctx);
}

namespace ExecuteObservable {
    export let PRINT_ERRORS_TO_STD_ERR = false;
}

export default ExecuteObservable