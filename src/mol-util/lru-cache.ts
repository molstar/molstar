
/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Adapted from LiteMol.
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LinkedList } from '../mol-data/generic';

export { LRUCache };

interface LRUCache<T> {
    entries: LinkedList<LRUCache.Entry<T>>,
    capacity: number
}

namespace LRUCache {
    export interface Entry<T> {
        key: string,
        data: T
    }

    function entry<T>(key: string, data: T): Entry<T> {
        return { key, data };
    }

    export function create<T>(capacity: number): LRUCache<T> {
        return {
            entries: LinkedList<Entry<T>>(),
            capacity: Math.max(1, capacity)
        };
    }

    export function get<T>(cache: LRUCache<T>, key: string) {
        for (let e = cache.entries.first; e; e = e.next) {
            if (e.value.key === key) {
                cache.entries.remove(e);
                cache.entries.addLast(e.value);
                return e.value.data;
            }
        }
        return void 0;
    }

    export function set<T>(cache: LRUCache<T>, key: string, data: T): T | undefined {
        let removed: T | undefined = undefined;
        if (cache.entries.count >= cache.capacity) {
            const first = cache.entries.first!;
            removed = first.value.data;
            cache.entries.remove(first);
        }
        cache.entries.addLast(entry(key, data));
        return removed;
    }
}