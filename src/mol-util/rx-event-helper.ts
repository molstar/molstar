/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Subject } from 'rxjs';

export class RxEventHelper {
    private _eventList: Subject<any>[] = [];
    create<T>() {
        const s = new Subject<T>();
        this._eventList.push(s);
        return s;
    }
    dispose() {
        for (const e of this._eventList) e.complete();
    }
}