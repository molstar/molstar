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
    private abortRequests = new Map<number, string | undefined>();
    private currentContext: { ctx: RuntimeContext, refCount: number }[] = [];

    readonly events = {
        progress: this.ev<TaskManager.ProgressEvent>(),
        finished: this.ev<{ id: number }>()
    };

    private track(internalId: number, taskId: number) {
        return (progress: Progress) => {
            if (progress.canAbort && progress.requestAbort && this.abortRequests.has(progress.root.progress.taskId)) {
                progress.requestAbort(this.abortRequests.get(taskId));
            }
            const elapsed = now() - progress.root.progress.startedTime;
            this.events.progress.next({
                id: internalId,
                level: elapsed < 250 ? 'none' : elapsed < 1500 ? 'background' : 'overlay',
                progress
            });
        };
    }

    async run<T>(task: Task<T>, createNewContext = false): Promise<T> {
        const id = this.id++;

        let ctx: TaskManager['currentContext'][0];

        if (createNewContext || this.currentContext.length === 0) {
            ctx = { ctx: CreateObservableCtx(task, this.track(id, task.id), 100), refCount: 1 };
        } else {
            ctx = this.currentContext[this.currentContext.length - 1];
            ctx.refCount++;
        }

        try {
            const ret = await ExecuteInContext(ctx.ctx, task);
            return ret;
        } finally {
            this.events.finished.next({ id });
            this.abortRequests.delete(task.id);
            ctx.refCount--;
            if (ctx.refCount === 0) arrayRemoveInPlace(this.currentContext, ctx);
        }
    }

    requestAbort(progress: Progress, reason?: string) {
        this.abortRequests.set(progress.root.progress.taskId, reason);
    }

    dispose() {
        this.ev.dispose();
    }
}

namespace TaskManager {
    export type ReportLevel = 'none' | 'background' | 'overlay'

    export interface ProgressEvent {
        id: number,
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