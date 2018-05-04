/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol
 * Copyright (c) 2016 - now David Sehnal, licensed under Apache 2.0, See LICENSE file for more info.
 */

import { Subject } from 'rxjs';
import { filter } from 'rxjs/operators';
import { Event } from '../event/event'

export class Dispatcher {
    LOG_DISPATCH_STREAM = false;

    private lanes: Subject<Event<any>>[] = [];
    constructor() {
        for (let i = 0; i <= Dispatcher.Lane.Job; i++) {
            this.lanes.push(new Subject<Event<any>>());
        }
    }

    dispatch<T>(event: Event<T>) {
        if (this.LOG_DISPATCH_STREAM) console.log(event.type.name, Dispatcher.Lane[event.type.lane], event.data);
        this.lanes[event.type.lane].next(event);
    }

    schedule(action: () => void, onError?: (e: string) => void, timeout = 1000 / 31) {
        return setTimeout(() => {
            if (onError) {
                try {
                    action.call(null)
                } catch (e) {
                    onError.call(null, '' + e);
                }
            } else {
                action.call(null);
            }
        }, timeout);
    }

    getStream<T>(type: Event.Type<T>): Event.Stream<T> {
        return this.lanes[type.lane].pipe(filter(e => e.type === type));
    }

    finished() {
        this.lanes.forEach(l => l.complete());
    }
}

export namespace Dispatcher {
    export enum Lane {
        Slow = 0,
        Fast = 1,
        Log = 2,
        Busy = 3,
        Transformer = 4,
        Job = 5
    }
}