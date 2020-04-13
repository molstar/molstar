/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { iterableToArray } from '../mol-data/util/array';

// TODO use set@@iterator when targeting es6

export namespace SetUtils {
    export function toArray<T>(set: ReadonlySet<T>) {
        return iterableToArray(set.values());
    }

    /** Test if set a contains all elements of set b. */
    export function isSuperset<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>) {
        let flag = true;
        setB.forEach(elem => {
            if (!setA.has(elem)) flag = false;
        });
        return flag;
    }

    /** Add all elements from `sets` to `out` */
    export function add<T>(out: Set<T>, ...sets: ReadonlySet<T>[]): Set<T> {
        for (let i = 0; i < sets.length; i++) {
            sets[i].forEach(elem => out.add(elem));
        }
        return out;
    }

    /** Create set containing elements of both set a and set b. */
    export function union<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        const union = new Set(setA);
        setB.forEach(elem => union.add(elem));
        return union;
    }

    export function unionMany<T>(...sets: ReadonlySet<T>[]) {
        if (sets.length === 0) return new Set<T>();
        if (sets.length === 1) new Set(sets[0]);
        const union = new Set(sets[0]);
        for (let i = 1, il = sets.length; i < il; i++) {
            sets[i].forEach(elem => union.add(elem));
        }
        return union;
    }

    export function unionManyArrays<T>(arrays: T[][]) {
        if (arrays.length === 0) return new Set<T>();
        const union = new Set(arrays[0]);
        for (let i = 1; i < arrays.length; i++) {
            for (const elem of arrays[i]) union.add(elem);
        }
        return union;
    }

    /** Create set containing elements of set a that are also in set b. */
    export function intersection<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        const intersection = new Set<T>();
        setB.forEach(elem => {
            if (setA.has(elem)) intersection.add(elem);
        });
        return intersection;
    }

    export function areIntersecting<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): boolean {
        let flag = false;
        setB.forEach(elem => {
            if (setA.has(elem)) flag = true;
        });
        return flag;
    }

    export function intersectionSize<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): number {
        let count = 0;
        setB.forEach(elem => {
            if (setA.has(elem)) count += 1;
        });
        return count;
    }

    /** Create set containing elements of set a that are not in set b. */
    export function difference<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        const difference = new Set(setA);
        setB.forEach(elem => difference.delete(elem));
        return difference;
    }

    /** Number of elements that are in set a but not in set b. */
    export function differenceSize<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): number {
        let count = setA.size;
        setA.forEach(elem => {
            if (setB.has(elem)) count -= 1;
        });
        return count;
    }

    /** Test if set a and b contain the same elements. */
    export function areEqual<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>) {
        if (setA.size !== setB.size) return false;
        let flag = true;
        setB.forEach(elem => {
            if (!setA.has(elem)) flag = false;
        });
        return flag;
    }
}