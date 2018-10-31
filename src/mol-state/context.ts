/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs'
import { StateObject } from './object';
import { Transform } from './transform';

interface StateContext {
    events: {
        object: {
            stateChanged: Subject<{ ref: Transform.Ref }>,
            propsChanged: Subject<{ ref: Transform.Ref, newProps: unknown }>,

            updated: Subject<{ ref: Transform.Ref, obj?: StateObject }>,
            replaced: Subject<{ ref: Transform.Ref, oldObj?: StateObject, newObj?: StateObject }>,
            created: Subject<{ ref: Transform.Ref, obj: StateObject }>,
            removed: Subject<{ ref: Transform.Ref, obj?: StateObject }>,
        },
        warn: Subject<string>
    },
    globalContext: unknown,
    defaultObjectProps: unknown
}

namespace StateContext {
    export function create(params: { globalContext: unknown, defaultObjectProps: unknown }): StateContext {
        return {
            events: {
                object: {
                    stateChanged: new Subject(),
                    propsChanged: new Subject(),
                    updated: new Subject(),
                    replaced: new Subject(),
                    created: new Subject(),
                    removed: new Subject()
                },
                warn: new Subject()
            },
            globalContext: params.globalContext,
            defaultObjectProps: params.defaultObjectProps
        }
    }
}

export { StateContext }