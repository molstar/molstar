/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export interface DefaultMap<K, V> extends Map<K, V> {
    /** Return the value for `key` when available or the default value otherwise. */
    getDefault: (key: K) => V
}

/** A `Map` instance with a `getDefault` method added. */
export function DefaultMap<K, V>(valueCtor: () => V) {
    const map = new Map<K, V>() as DefaultMap<K, V>;
    map.getDefault = (key: K) => {
        if (map.has(key)) return map.get(key)!;
        const value = valueCtor();
        map.set(key, value);
        return value;
    };
    return map;
}

// TODO currently not working, see https://github.com/Microsoft/TypeScript/issues/10853
// /** A `Map` with a `getDefault` method added. */
// export class DefaultMap<K, V> extends Map<K, V> {
//     constructor(private valueCtor: () => V, entries?: ReadonlyArray<[K, V]>) {
//         super(entries)
//     }
//     /** Return the value for `key` when available or the default value otherwise. */
//     getDefault(key: K) {
//         if (this.has(key)) return this.get(key)!
//         const value = this.valueCtor()
//         this.set(key, value)
//         return value
//     }
// }

export function arrayMapAdd<K, V>(map: Map<K, V[]>, key: K, value: V) {
    const a = map.get(key);
    if (a) {
        a.push(value);
    } else {
        map.set(key, [value]);
    }
}