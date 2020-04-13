/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { iterableToArray } from '../util';

// TODO: rename to "linear map" and just do key value mapping from index?

/** Immutable by convention IntMap */
interface IntMap<T> {
    has(key: number): boolean,
    keys(): IterableIterator<number>,
    values(): IterableIterator<T>,
    get(key: number): T,
    readonly size: number
}

namespace IntMap {
    export const Empty: IntMap<any> = new Map<number, any>();

    export interface Mutable<T> extends IntMap<T> {
        set(key: number, value: T): void;
    }

    export function keyArray<T>(map: IntMap<T>): number[] {
        return iterableToArray(map.keys());
    }

    export function Mutable<T>(): Mutable<T> {
        return new Map<number, T>() as Mutable<T>;
    }

    export function asImmutable<T>(map: IntMap<T>): IntMap<T> {
        return map;
    }

    export function copy<T>(map: IntMap<T>): Mutable<T> {
        const ret = Mutable<T>();
        const it = map.keys();
        while (true) {
            const { done, value } = it.next();
            if (done) break;
            ret.set(value, map.get(value));
        }
        return ret;
    }

    export function addFrom<T>(map: Mutable<T>, src: IntMap<T>) {
        const it = src.keys();
        while (true) {
            const { done, value } = it.next();
            if (done) break;
            map.set(value, src.get(value));
        }
        return map;
    }
}

export default IntMap;