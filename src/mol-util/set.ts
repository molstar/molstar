/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// TODO remove Array.from workaround when targeting ES6

export namespace SetUtils {
    /** Test if set a contains all elements of set b. */
    export function isSuperset<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>) {
        for (const elm of Array.from(setB)) {
            if (!setA.has(elm)) return false;
        }
        return true;
    }

    /** Add all elements from `sets` to `out` */
    export function add<T>(out: Set<T>, ...sets: ReadonlySet<T>[]): Set<T> {
        for (let i = 0; i < sets.length; i++) {
            for (const elem of Array.from(sets[i])) out.add(elem);
        }
        return out;
    }

    /** Create set containing elements of both set a and set b. */
    export function union<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        const union = new Set(setA);
        for (const elem of Array.from(setB)) union.add(elem);
        return union;
    }

    export function unionMany<T>(...sets: ReadonlySet<T>[]) {
        if (sets.length === 0) return new Set<T>();
        if (sets.length === 1) new Set(sets[0]);
        const union = new Set(sets[0]);
        for (let i = 1; i < sets.length; i++) {
            for (const elem of Array.from(sets[i])) union.add(elem);
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
        for (const elem of Array.from(setB)) {
            if (setA.has(elem)) intersection.add(elem);
        }
        return intersection;
    }

    export function areIntersecting<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): boolean {
        for (const elem of Array.from(setB)) {
            if (setA.has(elem)) return true;
        }
        return false;
    }

    /** Create set containing elements of set a that are not in set b. */
    export function difference<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        const difference = new Set(setA);
        for (const elem of Array.from(setB)) difference.delete(elem);
        return difference;
    }

    /** Test if set a and b contain the same elements. */
    export function areEqual<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>) {
        if (setA.size !== setB.size) return false
        for (const elm of Array.from(setB)) {
            if (!setA.has(elm)) return false;
        }
        return true;
    }
}