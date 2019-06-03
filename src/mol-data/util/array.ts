/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * from https://github.com/dsehnal/CIFTools.js
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ArrayCtor } from '../../mol-util/type-helpers';

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

/** Fills the array so that array[0] = start and array[array.length - 1] = end */
export function createRangeArray(start: number, end: number, ctor?: ArrayCtor<number>) {
    const len = end - start + 1;
    const array = ctor ? new ctor(len) : new Int32Array(len);
    for (let i = 0; i < len; i++) {
        array[i] = i + start;
    }
    return array;
}

export function arrayPickIndices<T>(array: ArrayLike<T>, indices: ArrayLike<number>) {
    const ret = new (arrayGetCtor(array))(indices.length);
    for (let i = 0, _i = indices.length; i < _i; i++) {
        ret[i] = array[indices[i]];
    }
    return ret;
}

export function arrayGetCtor<T>(data: ArrayLike<T>): ArrayCtor<T> {
    const ret = (data as any).constructor;
    if (!ret) throw new Error('data does not define a constructor and it should');
    return ret;
}