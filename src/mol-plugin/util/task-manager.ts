/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task, Progress } from 'mol-task';
import { RxEventHelper } from 'mol-util/rx-event-helper';
import { now } from 'mol-util/now';

export { TaskManager }

class TaskManager {
    private ev = RxEventHelper.create();
    private id = 0;
    private abortRequests = new Map<number, string | undefined>();

    readonly events = {
        progress: this.ev<TaskManager.ProgressEvent>(),
        finished: this.ev<{ id: number }>()
    };

    private track(internalId: number, taskId: number) {
        return (progress: Progress) => {
            if (progress.canAbort && progress.requestAbort && this.abortRequests.has(taskId)) {
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

    async run<T>(task: Task<T>): Promise<T> {
        const id = this.id++;
        try {
            const ret = await task.run(this.track(id, task.id), 100);
            return ret;
        } finally {
            this.events.finished.next({ id });
            this.abortRequests.delete(task.id);
        }
    }

    requestAbort(task: Task<any> | number, reason?: string) {
        this.abortRequests.set(typeof task === 'number' ? task : task.id, reason);
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
        })
    }
}