/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs'
import { StateObject } from './object';
import { Task } from 'mol-task';
import { Transform } from './tree/transform';

interface StateContext {
    events: {
        object: {
            stateChanged: Subject<{ ref: Transform.Ref }>,
            updated: Subject<{ ref: Transform.Ref }>,
            replaced: Subject<{ ref: Transform.Ref, old?: StateObject }>,
            created: Subject<{ ref: Transform.Ref }>,
            removed: Subject<{ ref: Transform.Ref }>,
        },
        warn: Subject<string>
    },
    globalContext: unknown,
    runTask<T>(task: T | Task<T>): T | Promise<T>
}

namespace StateContext {
    export function create(globalContext?: unknown/* task?: { observer?: Progress.Observer, updateRateMs?: number } */): StateContext {
        return {
            events: {
                object: {
                    stateChanged: new Subject(),
                    updated: new Subject(),
                    replaced: new Subject(),
                    created: new Subject(),
                    removed: new Subject()
                },
                warn: new Subject()
            },
            globalContext,
            runTask<T>(t: T | Task<T>) {
                if (typeof (t as any).run === 'function') return (t as Task<T>).run();
                return t as T;
            }
        }
    }
}

export { StateContext }