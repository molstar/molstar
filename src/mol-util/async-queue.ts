/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { arrayRemoveInPlace } from './array';
import { Subject } from 'rxjs';

export class AsyncQueue<T> {
    private queue: T[] = [];
    private signal = new Subject<{ v: T, stillPresent: boolean }>();

    get length() { return this.queue.length; }

    enqueue(v: T) {
        this.queue.push(v);
        if (this.queue.length === 1) return true;
        return this.waitFor(v);
    }

    handled(v: T) {
        arrayRemoveInPlace(this.queue, v);
        if (this.queue.length > 0) {
            this.signal.next({ v: this.queue[0], stillPresent: true });
        }
    }

    remove(v: T) {
        const rem = arrayRemoveInPlace(this.queue, v);
        if (rem)
            this.signal.next({ v, stillPresent: false });
        return rem;
    }

    private waitFor(t: T): Promise<boolean> {
        return new Promise(res => {
            const sub = this.signal.subscribe(({ v, stillPresent: removed }) => {
                if (v === t) {
                    sub.unsubscribe();
                    res(removed);
                }
            });
        });
    }
}