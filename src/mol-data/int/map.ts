/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { iterableToArray } from '../util'

/** Immutable by convention IntMap */
interface IntMap<T> {
    has(key: number): boolean,
    keys(): IterableIterator<number>,
    values(): IterableIterator<T>,
    get(key: number): T | undefined,
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

    export function ofObj<T>(obj: { [key: number]: T }): IntMap<T> {
        const keys = Object.keys(obj);
        const ret = new Map<number, T>();
        for (let i = 0, _i = keys.length; i < _i; i++) {
            const k = keys[i];
            ret.set(+k, obj[k as any]);
        }
        return ret as IntMap<T>;
    }

    export function Mutable<T>(): Mutable<T> {
        return new Map<number, T>() as Mutable<T>;
    }

    export function asImmutable<T>(map: IntMap<T>): IntMap<T> {
        return map;
    }
}

export default IntMap;