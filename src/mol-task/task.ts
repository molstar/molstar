/**
 * Copyright (c) 2017-2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeContext } from './execution/runtime-context';
import { Progress } from './execution/progress';
import { ExecuteObservable, ExecuteObservableChild, ExecuteInContext } from './execution/observable';
import { SyncRuntimeContext } from './execution/synchronous';
import { idFactory } from '../mol-util/id-factory';

/** A "named function wrapper" with built in "computation tree progress tracking". */
interface Task<T> {
    /** run the task without observation */
    run(): Promise<T>,
    /** run the task with the specified observer, default updateRate is 250ms */
    run(observer: Progress.Observer, updateRateMs?: number): Promise<T>,

    /**
     * Run a child task that adds a new node to the progress tree. Allows to passing the progress so
     * that the progress tree can be kept in a "good state" without having to separately call update.
     */
    runAsChild(ctx: RuntimeContext, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T>

    /** Run the task on the specified context. */
    runInContext(ctx: RuntimeContext): Promise<T>

    readonly id: number,
    readonly name: string
}

namespace Task {
    class Impl<T> implements Task<T> {
        readonly id: number;

        run(observer?: Progress.Observer, updateRateMs = 250): Promise<T> {
            if (observer) return ExecuteObservable(this, observer, updateRateMs);
            return this.f(SyncRuntimeContext);
        }

        runAsChild(ctx: RuntimeContext, progress?: string | Partial<RuntimeContext.ProgressUpdate>): Promise<T> {
            if (ctx.isSynchronous) return this.f(SyncRuntimeContext);
            return ExecuteObservableChild(ctx, this, progress as string | Partial<RuntimeContext.ProgressUpdate>);
        }

        runInContext(ctx: RuntimeContext): Promise<T> {
            if (ctx.isSynchronous) return this.f(SyncRuntimeContext);
            return ExecuteInContext(ctx, this);
        }

        constructor(public name: string, public f: (ctx: RuntimeContext) => Promise<T>, public onAbort?: () => void) {
            this.id = getNextId();
        }
    }

    export function is<T = any>(t: any): t is Task<T> {
        const _t = t as Task<any>;
        return !!t && typeof _t.id === 'number' && typeof _t.name === 'string' && !!_t.run;
    }

    export interface Aborted { isAborted: true, reason: string, toString(): string }
    export function isAbort(e: any): e is Aborted { return !!e && !!e.isAborted; }
    export function Aborted(reason: string): Aborted { return { isAborted: true, reason, toString() { return `Aborted${reason ? ': ' + reason : ''}`; } }; }

    export function create<T>(name: string, f: (ctx: RuntimeContext) => Promise<T>, onAbort?: () => void): Task<T> {
        return new Impl(name, f, onAbort);
    }

    export function constant<T>(name: string, value: T): Task<T> { return create(name, async ctx => value); }
    export function empty(): Task<void> { return create('', async ctx => {}); }
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

    const getNextId = idFactory(0, 0x3fffffff);
}

export { Task };