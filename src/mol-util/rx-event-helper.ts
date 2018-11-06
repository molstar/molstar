/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs';

export { RxEventHelper }

interface RxEventHelper {
    <T>(): Subject<T>,
    dispose(): void
}

namespace RxEventHelper {
    export function create(): RxEventHelper {
        const helper = new _RxEventHelper();
        const ret: RxEventHelper = (<T>() => helper.create<T>()) as RxEventHelper;
        ret.dispose = () => helper.dispose();
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
    dispose() {
        if (this._disposed) return;
        for (const e of this._eventList) e.complete();
        this._disposed = true;
    }
}