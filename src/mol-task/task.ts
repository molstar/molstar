/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RuntimeContext } from './execution/runtime-context'

// A "named function wrapper" with built in "computation tree progress tracking".
// Use Run(t, ?observer, ?updateRate) to execute
interface Task<T> {
    readonly id: number,
    readonly name: string,
    // Do not call this directly, use Run.
    readonly __f: (ctx: RuntimeContext) => Promise<T>,
    // Do not call this directly, use Run.
    readonly __onAbort: (() => void) | undefined
}

namespace Task {
    export interface Aborted { isAborted: true, reason: string }
    export function isAbort(e: any): e is Aborted { return !!e && !!e.isAborted; }
    export function Aborted(reason: string): Aborted { return { isAborted: true, reason }; }

    export function create<T>(name: string, f: (ctx: RuntimeContext) => Promise<T>, onAbort?: () => void): Task<T> {
        return { id: nextId(), name, __f: f, __onAbort: onAbort };
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