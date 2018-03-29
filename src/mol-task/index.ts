/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from './task'
import { RuntimeContext } from './execution/runtime-context'
import { ExecuteSynchronous } from './execution/synchronous'
import { ExecuteObservable } from './execution/observable'
import { Progress } from './execution/progress'
import { now } from './util/now'
import { Scheduler } from './util/scheduler'
import { MultistepTask } from './util/multistep'
import { ChunkedSubtask } from './util/chunked'

// Run the task without the ability to observe its progress.
function Run<T>(task: Task<T>): Promise<T>;
// Run the task with the ability to observe its progress and request cancellation.
function Run<T>(task: Task<T>, observer: Progress.Observer, updateRateMs?: number): Promise<T>;
function Run<T>(task: Task<T>, observer?: Progress.Observer, updateRateMs?: number): Promise<T> {
    if (observer) return ExecuteObservable(task, observer, updateRateMs || 250);
    return ExecuteSynchronous(task);
}

export { Task, RuntimeContext, Progress, Run, now, Scheduler, MultistepTask, ChunkedSubtask }