/**
 * Copyright (c) 2018-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// make use of https://github.com/tc39/proposal-set-methods once available

export namespace SetUtils {
    export function toArray<T>(set: ReadonlySet<T>) {
        return Array.from(set.values());
    }

    /** Test if set a contains all elements of set b. */
    export function isSuperset<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>) {
        if (setA === setB) return true;
        if (setA.size < setB.size) return false;
        for (const elem of setB) {
            if (!setA.has(elem)) return false;
        }
        return true;
    }

    /** Add all elements from `sets` to `out` */
    export function add<T>(out: Set<T>, ...sets: ReadonlySet<T>[]): Set<T> {
        for (let i = 0; i < sets.length; i++) {
            for (const elem of sets[i]) out.add(elem);
        }
        return out;
    }

    /** Create set containing elements of both set a and set b. */
    export function union<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        const union = new Set(setA);
        if (setA === setB) return union;
        for (const elem of setB) union.add(elem);
        return union;
    }

    export function unionMany<T>(...sets: ReadonlySet<T>[]) {
        if (sets.length === 0) return new Set<T>();
        const union = new Set(sets[0]);
        for (let i = 1, il = sets.length; i < il; i++) {
            for (const elem of sets[i]) union.add(elem);
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
        if (setA === setB) return new Set(setA);
        const intersection = new Set<T>();
        for (const elem of setB) {
            if (setA.has(elem)) intersection.add(elem);
        }
        return intersection;
    }

    export function areIntersecting<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): boolean {
        if (setA === setB) return setA.size > 0;
        if (setA.size < setB.size) [setA, setB] = [setB, setA];
        for (const elem of setB) {
            if (setA.has(elem)) return true;
        }
        return false;
    }

    export function intersectionSize<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): number {
        if (setA === setB) return setA.size;
        if (setA.size < setB.size) [setA, setB] = [setB, setA];
        let count = 0;
        for (const elem of setB) {
            if (setA.has(elem)) count += 1;
        }
        return count;
    }

    /** Create set containing elements of set a that are not in set b. */
    export function difference<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): Set<T> {
        if (setA === setB) return new Set();
        const difference = new Set(setA);
        for (const elem of setB) difference.delete(elem);
        return difference;
    }

    /** Number of elements that are in set a but not in set b. */
    export function differenceSize<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>): number {
        if (setA === setB) return 0;
        let count = setA.size;
        for (const elem of setA) {
            if (setB.has(elem)) count -= 1;
        }
        return count;
    }

    /** Test if set a and b contain the same elements. */
    export function areEqual<T>(setA: ReadonlySet<T>, setB: ReadonlySet<T>) {
        if (setA === setB) return true;
        if (setA.size !== setB.size) return false;
        for (const elem of setB) {
            if (!setA.has(elem)) return false;
        }
        return true;
    }
}