/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeContext } from './execution/runtime-context'
import { Progress } from './execution/progress'
import { ExecuteObservable, ExecuteObservableChild } from './execution/observable';
import { SyncRuntimeContext } from 'mol-task/execution/synchronous';

// A "named function wrapper" with built in "computation tree progress tracking".
// Use Run(t, ?observer, ?updateRate) to execute
interface Task<T> {
    // run the task without observation
    run(): Promise<T>,
    // run the task with the specified observer, default updateRate is 250ms
    run(observer: Progress.Observer, updateRateMs?: number): Promise<T>,

    // Run a child task that adds a new node to the progress tree.
    // Allow to pass the progress so that the progress tree can be kept in a "good state" without having to separately call update.
    runAsChild(ctx: RuntimeContext, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T>

    readonly id: number,
    readonly name: string
}

namespace Task {
    class Impl<T> implements Task<T> {
        readonly id: number;

        run(observer?: Progress.Observer, updateRateMs?: number): Promise<T> {
            if (observer) return ExecuteObservable(this, observer, updateRateMs as number || 250);
            return this.f(SyncRuntimeContext);
        }

        runAsChild(ctx: RuntimeContext, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T> {
            if (ctx.isSynchronous) return this.f(SyncRuntimeContext);
            return ExecuteObservableChild(ctx, this, progress as string | Partial<RuntimeContext.ProgressUpdate>);
        }

        constructor(public name: string, public f: (ctx: RuntimeContext) => Promise<T>, public onAbort?: () => void) {
            this.id = nextId();
        }
    }

    export interface Aborted { isAborted: true, reason: string }
    export function isAbort(e: any): e is Aborted { return !!e && !!e.isAborted; }
    export function Aborted(reason: string): Aborted { return { isAborted: true, reason }; }

    export function create<T>(name: string, f: (ctx: RuntimeContext) => Promise<T>, onAbort?: () => void): Task<T> {
        return new Impl(name, f, onAbort);
    }

    export function constant<T>(name: string, value: T): Task<T> { return create(name, async ctx => value); }
    export function fail(name: string, reason: string): Task<any> { return create(name, async ctx => { throw new Error(reason); }); }

    export interface Progress {
        taskId: number,
        taskName: string,
        startedTime: number,
        message: string,
        canAbort: boolean,
        isIndeterminate: boolean,
        current: number,
        max: number
    }

    let _id = 0;
    function nextId() {
        const ret = _id;
        _id = (_id + 1) % 0x3fffffff;
        return ret;
    }
}

export { Task }