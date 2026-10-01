/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject, BehaviorSubject } from 'rxjs';

export { RxEventHelper };

interface RxEventHelper {
    <T>(): Subject<T>,
    behavior<T>(v: T): BehaviorSubject<T>,
    dispose(): void
}

namespace RxEventHelper {
    export function create(): RxEventHelper {
        const helper = new _RxEventHelper();
        const ret: RxEventHelper = (<T>() => helper.create<T>()) as RxEventHelper;
        ret.dispose = () => helper.dispose();
        ret.behavior = (v) => helper.behavior(v);
        return ret;
    }
}

class _RxEventHelper {
    private _eventList: Subject<any>[] = [];
    private _disposed = false;

    create<T>() {
        const s = new Subject<T>();
        this._eventList.push(s);
        return s;
    }

    behavior<T>(v: T) {
        const s = new BehaviorSubject<T>(v);
        this._eventList.push(s);
        return s;
    }

    dispose() {
        if (this._disposed) return;
        for (const e of this._eventList) e.complete();
        this._disposed = true;
    }
}