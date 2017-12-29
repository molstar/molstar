/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

export function arrayFind<T>(array: ArrayLike<T>, f: (v: T) => boolean): T | undefined {
    for (let i = 0, _i = array.length; i < _i; i++) {
        if (f(array[i])) return array[i];
    }
    return void 0;
}

export function iterableToArray<T>(it: IterableIterator<T>): T[] {
    if (Array.from) return Array.from(it);

    const ret = [];
    while (true) {
        const { done, value } = it.next();
        if (done) break;
        ret[ret.length] = value;
    }
    return ret;
}