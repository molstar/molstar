/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Observable } from 'rxjs';
import { Context } from '../context/context'
import { Dispatcher } from '../service/dispatcher'

export interface Event<T> {
    type: Event.Type<T>;
    data: T;
}

export namespace Event {
    export type Stream<T> = Observable<Event<T>>;

    import Lane = Dispatcher.Lane

    export type Any = Event<any>
    export type AnyType = Type<any>

    export interface Type<T> {
        name: string,
        lane: Lane,
        dispatch(context: Context, data: T): void;
        getStream(context: Context): Stream<T>;
    }

    const EventPrototype = {
        dispatch<T>(this: any, context: Context, data: T) { context.dispatcher.dispatch({ type: this, data }) },
        getStream(this: any, context: Context) { return context.dispatcher.getStream(this); }
    }

    export function create<T>(name: string, lane: Dispatcher.Lane): Type<T> {
        return Object.create(EventPrototype, {
            name: { writable: false, configurable: false, value: name },
            lane: { writable: false, configurable: false, value: lane }
        });
    }
}