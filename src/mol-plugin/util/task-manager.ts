/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, Progress, RuntimeContext } from '../../mol-task';
import { RxEventHelper } from '../../mol-util/rx-event-helper';
import { now } from '../../mol-util/now';
import { CreateObservableCtx, ExecuteInContext } from '../../mol-task/execution/observable';
import { arrayRemoveInPlace } from '../../mol-util/array';

export { TaskManager };

class TaskManager {
    private ev = RxEventHelper.create();
    private id = 0;
    private runningTasks = new Set<number>();
    private abortRequests = new Map<number, string | undefined>();
    private options = new Map<number, { useOverlay: boolean }>();
    private currentContext: { ctx: RuntimeContext, refCount: number }[] = [];

    readonly events = {
        progress: this.ev<TaskManager.ProgressEvent>(),
        finished: this.ev<{ id: number }>()
    };

    private tryGetAbortTaskId(node: Progress.Node): number | undefined {
        if (this.abortRequests.has(node.progress.taskId)) return node.progress.taskId;
        for (const c of node.children) {
            const abort = this.tryGetAbortTaskId(c);
            if (abort !== void 0) return abort;
        }
        return void 0;
    }

    private track(internalId: number, taskId: number) {
        return (progress: Progress) => {
            if (progress.canAbort && progress.requestAbort) {
                const abortTaskId = this.tryGetAbortTaskId(progress.root);
                if (abortTaskId !== void 0) progress.requestAbort(this.abortRequests.get(abortTaskId));
            }
            const elapsed = now() - progress.root.progress.startedTime;
            this.events.progress.next({
                id: internalId,
                useOverlay: this.options.get(taskId)?.useOverlay,
                level: elapsed < 250 ? 'none' : 'background',
                progress
            });
        };
    }

    async run<T>(task: Task<T>, params?: { createNewContext?: boolean, useOverlay?: boolean }): Promise<T> {
        const id = this.id++;

        let ctx: TaskManager['currentContext'][0];

        if (params?.createNewContext || this.currentContext.length === 0) {
            ctx = { ctx: CreateObservableCtx(task, this.track(id, task.id), 100), refCount: 1 };
        } else {
            ctx = this.currentContext[this.currentContext.length - 1];
            ctx.refCount++;
        }

        try {
            this.options.set(task.id, { useOverlay: !!params?.useOverlay });
            this.runningTasks.add(task.id);
            const ret = await ExecuteInContext(ctx.ctx, task);
            return ret;
        } finally {
            this.options.delete(task.id);
            this.runningTasks.delete(task.id);
            this.events.finished.next({ id });
            this.abortRequests.delete(task.id);
            ctx.refCount--;
            if (ctx.refCount === 0) arrayRemoveInPlace(this.currentContext, ctx);
        }
    }

    requestAbortAll(reason?: string) {
        this.runningTasks.forEach(id => this.abortRequests.set(id, reason));
    }

    requestAbort(taskIdOrProgress: number | Progress, reason?: string) {
        const id = typeof taskIdOrProgress === 'number'
            ? taskIdOrProgress
            : taskIdOrProgress.root.progress.taskId;
        this.abortRequests.set(id, reason);
    }

    dispose() {
        this.ev.dispose();
    }
}

namespace TaskManager {
    export type ReportLevel = 'none' | 'background'

    export interface ProgressEvent {
        id: number,
        useOverlay?: boolean,
        level: ReportLevel,
        progress: Progress
    }

    function delay(time: number): Promise<void> {
        return new Promise(res => setTimeout(res, time));
    }
    export function testTask(N: number) {
        return Task.create('Test', async ctx => {
            let i = 0;
            while (i < N) {
                await delay(100 + Math.random() * 200);
                if (ctx.shouldUpdate) {
                    await ctx.update({ message: 'Step ' + i, current: i, max: N, isIndeterminate: false });
                }
                i++;
            }
        });
    }
}