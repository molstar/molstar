/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

interface UniqueArray<K, T = K> {
    keys: Set<K>,
    array: T[]
}

namespace UniqueArray {
    export function create<K, T = K>(): UniqueArray<K, T> {
        return { keys: new Set<K>(), array: [] };
    }

    export function add<K, T>({ keys, array }: UniqueArray<K, T>, key: K, value: T) {
        if (keys.has(key)) return false;
        keys.add(key);
        array[array.length] = value;
        return true;
    }

    export function has<K, T>({ keys }: UniqueArray<K, T>, key: K) {
        return keys.has(key);
    }
}

export { UniqueArray };