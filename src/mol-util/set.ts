/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// TODO remove Array.from workaround when targeting ES6

/** Test if set a contains all elements of set b. */
export function isSuperset<T>(setA: Set<T>, setB: Set<T>) {
    for (const elm of Array.from(setB)) {
        if (!setA.has(elm)) return false;
    }
    return true;
}

/** Create set containing elements of both set a and set b. */
export function union<T>(setA: Set<T>, setB: Set<T>): Set<T> {
    const union = new Set(setA);
    for (const elem of Array.from(setB)) union.add(elem);
    return union;
}

/** Create set containing elements of set a that are also in set b. */
export function intersection<T>(setA: Set<T>, setB: Set<T>): Set<T> {
    const intersection = new Set();
    for (const elem of Array.from(setB)) {
        if (setA.has(elem)) intersection.add(elem);
    }
    return intersection;
}

/** Create set containing elements of set a that are not in set b. */
export function difference<T>(setA: Set<T>, setB: Set<T>): Set<T> {
    const difference = new Set(setA);
    for (const elem of Array.from(setB)) difference.delete(elem);
    return difference;
}