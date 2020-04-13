/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from '../task';
import { RuntimeContext } from './runtime-context';
import { Progress } from './progress';
import { now } from '../../mol-util/now';
import { Scheduler } from '../util/scheduler';
import { UserTiming } from '../util/user-timing';
import { isDebugMode } from '../../mol-util/debug';

interface ExposedTask<T> extends Task<T> {
    f: (ctx: RuntimeContext) => Promise<T>,
    onAbort?: () => void
}

export function ExecuteObservable<T>(task: Task<T>, observer: Progress.Observer, updateRateMs = 250) {
    const info = ProgressInfo(task, observer, updateRateMs);
    const ctx = new ObservableRuntimeContext(info, info.root);
    return execute(task as ExposedTask<T>, ctx);
}

export function CreateObservableCtx<T>(task: Task<T>, observer: Progress.Observer, updateRateMs = 250) {
    const info = ProgressInfo(task, observer, updateRateMs);
    return new ObservableRuntimeContext(info, info.root);
}

export function ExecuteInContext<T>(ctx: RuntimeContext, task: Task<T>) {
    return execute(task as ExposedTask<T>, ctx as ObservableRuntimeContext);
}

export function ExecuteObservableChild<T>(ctx: RuntimeContext, task: Task<T>, progress?: string | Partial<RuntimeContext.ProgressUpdate>) {
    return (ctx as ObservableRuntimeContext).runChild(task as ExposedTask<T>, progress);
}

function defaultProgress(task: Task<any>): Task.Progress {
    return {
        taskId: task.id,
        taskName: task.name,
        message: '',
        startedTime: 0,
        canAbort: true,
        isIndeterminate: true,
        current: 0,
        max: 0
    };
}

interface ProgressInfo {
    updateRateMs: number,
    lastNotified: number,
    observer: Progress.Observer,

    abortToken: { abortRequested: boolean, treeAborted: boolean, reason: string },
    taskId: number;
    root: Progress.Node;
    tryAbort: (reason?: string) => void;
}

function ProgressInfo(task: Task<any>, observer: Progress.Observer, updateRateMs: number): ProgressInfo {
    const abortToken: ProgressInfo['abortToken'] = { abortRequested: false, treeAborted: false, reason: '' };

    return {
        updateRateMs,
        lastNotified: now(),
        observer,
        abortToken,
        taskId: task.id,
        root: { progress: defaultProgress(task), children: [] },
        tryAbort: createAbortFunction(abortToken)
    };
}

function createAbortFunction(token: ProgressInfo['abortToken']) {
    return (reason?: string) => {
        token.abortRequested = true;
        token.reason = reason || token.reason;
    };
}

function cloneTree(root: Progress.Node): Progress.Node {
    return { progress: { ...root.progress }, children: root.children.map(cloneTree) };
}

function canAbort(root: Progress.Node): boolean {
    return root.progress.canAbort && root.children.every(canAbort);
}

function snapshotProgress(info: ProgressInfo): Progress {
    return { root: cloneTree(info.root), canAbort: canAbort(info.root), requestAbort: info.tryAbort };
}

async function execute<T>(task: ExposedTask<T>, ctx: ObservableRuntimeContext) {
    UserTiming.markStart(task);
    ctx.node.progress.startedTime = now();
    try {
        const ret = await task.f(ctx);
        UserTiming.markEnd(task);
        UserTiming.measure(task);
        if (ctx.info.abortToken.abortRequested) {
            abort(ctx.info);
        }
        return ret;
    } catch (e) {
        if (Task.isAbort(e)) {
            ctx.isAborted = true;

            // wait for all child computations to go thru the abort phase.
            if (ctx.node.children.length > 0) {
                await new Promise(res => { ctx.onChildrenFinished = res; });
            }
            if (task.onAbort) {
                task.onAbort();
            }
        }
        if (isDebugMode) console.error(e);
        throw e;
    }
}

function abort(info: ProgressInfo) {
    if (!info.abortToken.treeAborted) {
        info.abortToken.treeAborted = true;
        abortTree(info.root);
        notifyObserver(info, now());
    }
    throw Task.Aborted(info.abortToken.reason);
}

function abortTree(root: Progress.Node) {
    const progress = root.progress;
    progress.isIndeterminate = true;
    progress.canAbort = false;
    progress.message = 'Aborting...';
    for (const c of root.children) abortTree(c);
}

// function shouldNotify(info: ProgressInfo, time: number) {
//     return time - info.lastNotified > info.updateRateMs;
// }

function notifyObserver(info: ProgressInfo, time: number) {
    info.lastNotified = time;
    const snapshot = snapshotProgress(info);
    info.observer(snapshot);
}

class ObservableRuntimeContext implements RuntimeContext {
    isSynchronous = false;

    isExecuting = true;
    lastUpdatedTime = 0;

    isAborted?: boolean;

    node: Progress.Node;
    info: ProgressInfo;

    // used for waiting for cancelled computation trees
    onChildrenFinished?: () => void = void 0;

    private checkAborted() {
        if (this.info.abortToken.abortRequested) {
            this.isAborted = true;
            abort(this.info);
        }
    }

    get shouldUpdate(): boolean {
        this.checkAborted();
        return now() - this.lastUpdatedTime > this.info.updateRateMs;
    }

    private updateProgress(update?: string | Partial<RuntimeContext.ProgressUpdate>) {
        this.checkAborted();

        if (!update) return;

        const progress = this.node.progress;
        if (typeof update === 'string') {
            progress.message = update;
            progress.isIndeterminate = true;
        } else {
            if (typeof update.canAbort !== 'undefined') progress.canAbort = update.canAbort;
            if (typeof update.message !== 'undefined') progress.message = update.message;
            if (typeof update.current !== 'undefined') progress.current = update.current;
            if (typeof update.max !== 'undefined') progress.max = update.max;
            progress.isIndeterminate = typeof progress.current === 'undefined' || typeof progress.max === 'undefined';
            if (typeof update.isIndeterminate !== 'undefined') progress.isIndeterminate = update.isIndeterminate;
        }
    }

    update(progress?: string | Partial<RuntimeContext.ProgressUpdate>, dontNotify?: boolean): Promise<void> | void {
        // The progress tracking and observer notification are separated
        // because the computation can have a tree structure.
        // All nodes of the tree should be regualarly updated at the specified frequency,
        // however, the notification should only be invoked once per the whole tree.

        this.lastUpdatedTime = now();
        this.updateProgress(progress);

        // TODO: do the shouldNotify check here?
        if (!!dontNotify /* || !shouldNotify(this.info, this.lastUpdatedTime)*/) return;

        notifyObserver(this.info, this.lastUpdatedTime);

        // The computation could have been aborted during the notifycation phase.
        this.checkAborted();

        return Scheduler.immediatePromise();
    }

    async runChild<T>(task: ExposedTask<T>, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T> {
        this.updateProgress(progress);

        // Create a new child context and add it to the progress tree.
        // When the child task finishes, remove the tree node.

        const node: Progress.Node = { progress: defaultProgress(task), children: [] };
        const children = this.node.children as Progress.Node[];
        children.push(node);
        const ctx = new ObservableRuntimeContext(this.info, node);
        try {
            return await execute(task as ExposedTask<T>, ctx);
        } catch (e) {
            if (Task.isAbort(e)) {
                // need to catch the error here because otherwise
                // promises for running child tasks in a tree-like computation
                // will get orphaned and cause "uncaught error in Promise".
                if (this.isAborted) return void 0 as any;
            }
            throw e;
        } finally {
            // remove the progress node after the computation has finished.
            const idx = children.indexOf(node);
            if (idx >= 0) {
                for (let i = idx, _i = children.length - 1; i < _i; i++) {
                    children[i] = children[i + 1];
                }
                children.pop();
            }
            if (children.length === 0 && this.onChildrenFinished) this.onChildrenFinished();
        }
    }

    constructor(info: ProgressInfo, node: Progress.Node) {
        this.node = node;
        this.info = info;
    }
}