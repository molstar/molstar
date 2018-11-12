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

    readonly events = {
        progress: this.ev<TaskManager.ProgressEvent>(),
        finished: this.ev<{ id: number }>()
    };

    private track(id: number) {
        return (progress: Progress) => {
            const elapsed = now() - progress.root.progress.startedTime;
            progress.root.progress.startedTime
            this.events.progress.next({
                id,
                level: elapsed < 250 ? 'none' : elapsed < 1500 ? 'background' : 'overlay',
                progress
            });
        };
    }

    async run<T>(task: Task<T>): Promise<T> {
        const id = this.id++;
        try {
            const ret = await task.run(this.track(id), 100);
            return ret;
        } finally {
            this.events.finished.next({ id });
        }
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