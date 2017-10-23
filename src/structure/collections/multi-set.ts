/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from './ordered-set'
import Iterator from './iterator'
import IntPair from './int-pair'
import { sortArray } from './sort'
import { hash1 } from './hash-functions'

type MultiSetElements = { [id: number]: OrderedSet, size: number, hashCode: number, keys: OrderedSet }
type MultiSet = number | MultiSetElements

namespace MultiSet {
    export const Empty: MultiSet = { size: 0, hashCode: 0, keys: OrderedSet.Empty };

    export function create(data: number | ArrayLike<number> | IntPair | { [id: number]: OrderedSet }): MultiSet {
        if (typeof data === 'number') return data;
        if (IntPair.is(data)) return IntPair.pack(data);
        if (isArrayLike(data)) return ofPackedPairs(data);
        return ofObject(data);
    }

    export function keys(set: MultiSet): OrderedSet {
        if (typeof set === 'number') return OrderedSet.ofSingleton(set);
        return set.keys;
    }

    export function hasKey(set: MultiSet, key: number): boolean {
        if (typeof set === 'number') return IntPair.fst(set) === key;
        return set.keys.has(key);
    }

    export function get(set: MultiSet, key: number): OrderedSet {
        if (!hasKey(set, key)) return OrderedSet.Empty;
        if (typeof set === 'number') return OrderedSet.ofSingleton(IntPair.snd(set));
        return set[key];
    }

    export function size(set: MultiSet) {
        if (typeof set === 'number') return 0;
        return set.size;
    }

    export function hashCode(set: MultiSet) {
        if (typeof set === 'number') return IntPair.packedHashCode(set);
        if (set.hashCode !== -1) return set.hashCode;
        return computeHash(set);
    }

    export function areEqual(a: MultiSet, b: MultiSet): boolean {
        if (a === b) return true;
        if (typeof a === 'number') {
            if (typeof b === 'number') return a === b;
            return false;
        }
        if (typeof b === 'number') return false;
        return areEqualEE(a, b);
    }

    export function areIntersecting(a: MultiSet, b: MultiSet): boolean {
        if (a === b) return true;
        if (typeof a === 'number') {
            if (typeof b === 'number') return a === b;
            return areIntersectingNE(a, b);
        }
        if (typeof b === 'number') return areIntersectingNE(b, a);
        return areIntersectingEE(a, b);
    }

    export function intersect(a: MultiSet, b: MultiSet): MultiSet {
        if (a === b) return a;
        if (typeof a === 'number') {
            if (typeof b === 'number') return a === b ? a : Empty;
            return intersectNE(a, b);
        }
        if (typeof b === 'number') return intersectNE(b, a);
        return intersectEE(a, b);
    }

    export function subtract(a: MultiSet, b: MultiSet): MultiSet {
        if (a === b) return Empty;
        if (typeof a === 'number') {
            if (typeof b === 'number') return a === b ? Empty : a;
            return subtractNE(a, b);
        }
        if (typeof b === 'number') return subtractEN(a, b);
        return subtractEE(a, b);
    }

    // TODO: union

    export function union(sets: ArrayLike<MultiSet>): MultiSet {
        return 0 as any;
    }

    class ElementsIterator implements Iterator<IntPair> {
        private pair = IntPair.zero();
        private unit = 0;
        private currentIndex = -1;
        private currentIterator: Iterator<number> = Iterator.Empty;

        [Symbol.iterator]() { return new ElementsIterator(this.elements); };
        done: boolean;
        next() { const value = this.move(); return { value, done: this.done } }

        move() {
            if (this.done) return this.pair;

            let next = this.currentIterator.move();
            if (this.currentIterator.done) {
                if (!this.advanceIterator()) {
                    this.done = true;
                    return this.pair;
                }
                next = this.currentIterator.move();
            }

            this.pair.snd = next;
            return this.pair;
        }

        private advanceIterator() {
            const keys = this.elements.keys;
            if (++this.currentIndex >= keys.size) return false;
            this.unit = keys.elementAt(this.currentIndex);
            this.pair.fst = this.unit;
            this.currentIterator = this.elements[this.unit].elements();
            return true;
        }

        constructor(private elements: MultiSetElements) {
            this.done = elements.keys.size === 0;
            this.advanceIterator();
        }
    }

    export function values(set: MultiSet): Iterator<IntPair> {
        if (typeof set === 'number') return Iterator.Value(IntPair.unpack1(set));
        return new ElementsIterator(set);
    }
}

const pair = IntPair.zero();


function isArrayLike(x: any): x is ArrayLike<number> {
    return x && (typeof x.length === 'number' && (x instanceof Array || !!x.buffer));
}

function ofObject(data: { [id: number]: OrderedSet }) {
    const keys = [];
    for (const _k of Object.keys(data)) {
        const k = +_k;
        if (data[k].size > 0) keys[keys.length] = k;
    }
    if (!keys.length) return MultiSet.Empty;
    if (keys.length === 1) {
        const set = data[keys[0]];
        if (set.size === 1) return IntPair.pack1(keys[0], set.elementAt(0));
    }
    return ofObject1(keys, data);
}

function ofObject1(keys: number[], data: { [id: number]: OrderedSet }) {
    if (keys.length === 1) {
        const k = keys[0];
        const set = data[k];
        if (set.size === 1) return IntPair.pack1(k, set.elementAt(0));
    }
    sortArray(keys);
    return _createObjectOrdered(OrderedSet.ofSortedArray(keys), data);
}

function ofObjectOrdered(keys: OrderedSet, data: { [id: number]: OrderedSet }) {
    if (keys.size === 1) {
        const k = keys.elementAt(0);
        const set = data[k];
        if (set.size === 1) return IntPair.pack1(k, set.elementAt(0));
    }
    return _createObjectOrdered(keys, data);
}

function _createObjectOrdered(keys: OrderedSet, data: { [id: number]: OrderedSet }) {
    const ret: MultiSetElements = Object.create(null);
    ret.keys = keys;
    let size = 0;
    for (let i = 0, _i = keys.size; i < _i; i++) {
        const k = keys.elementAt(i);
        const set = data[k];
        ret[k] = set;
        size += set.size;
    }
    ret.size = size;
    ret.hashCode = -1;
    return ret;
}

function ofPackedPairs(xs: ArrayLike<number>): MultiSet {
    if (xs.length === 0) return MultiSet.Empty;
    const sets: { [key: number]: number[] } = Object.create(null);
    const p = IntPair.zero();
    for (let i = 0, _i = xs.length; i < _i; i++) {
        IntPair.unpack(xs[i], p);
        const set = sets[p.fst];
        if (set) set[set.length] = p.snd;
        else sets[p.fst] = [p.snd];
    }
    const ret: { [key: number]: OrderedSet } = Object.create(null);
    const keys = [];
    for (const _k of Object.keys(sets)) {
        const k = +_k;
        keys[keys.length] = k;
        ret[k] = OrderedSet.ofSortedArray(sortArray(sets[k]));
    }
    return ofObject1(keys, ret);
}

function computeHash(set: MultiSetElements) {
    const { keys } = set;
    let hash = 23;
    for (let i = 0, _i = keys.size; i < _i; i++) {
        const k = keys.elementAt(i);
        hash = (31 * hash + k) | 0;
        hash = (31 * hash + OrderedSet.hashCode(set[k])) | 0;
    }
    hash = (31 * hash + set.size) | 0;
    hash = hash1(hash);
    set.hashCode = hash;
    return hash;
}

function areEqualEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return true;
    if (a.size !== b.size) return false;

    const keys = a.keys;
    if (!OrderedSet.areEqual(keys, b.keys)) return false;
    for (let i = 0, _i = keys.size; i < _i; i++) {
        const k = keys.elementAt(i);
        if (!OrderedSet.areEqual(a[k], b[k])) return false;
    }
    return true;
}

function areIntersectingNE(a: number, b: MultiSetElements) {
    IntPair.unpack(a, pair);
    return b.keys.has(pair.fst) && b[pair.fst].has(pair.snd);
}

function areIntersectingEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return true;
    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return false;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));
    for (let i = start; i < end; i++) {
        const k = keysA.elementAt(i);
        if (keysB.has(k) && OrderedSet.areIntersecting(a[k], b[k])) return true;
    }
    return false;
}

function intersectNE(a: number, b: MultiSetElements) {
    IntPair.unpack(a, pair);
    return b.keys.has(pair.fst) && b[pair.fst].has(pair.snd) ? a : MultiSet.Empty;
}

function intersectEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return MultiSet.Empty;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));

    const keys = [], ret = Object.create(null);
    for (let i = start; i < end; i++) {
        const k = keysA.elementAt(i);
        if (keysB.has(k)) {
            const intersection = OrderedSet.intersect(a[k], b[k]);
            if (intersection.size > 0) {
                keys[keys.length] = k;
                ret[k] = intersection;
            }
        }
    }
    return ofObjectOrdered(OrderedSet.ofSortedArray(keys), ret);
}

function subtractNE(a: number, b: MultiSetElements) {
    IntPair.unpack(a, pair);
    return b.keys.has(pair.fst) && b[pair.fst].has(pair.snd) ? MultiSet.Empty : a;
}

function subtractEN(a: MultiSetElements, b: number): MultiSet {
    const aKeys =  a.keys;
    IntPair.unpack(b, pair);
    if (!aKeys.has(pair.fst) || !a[pair.fst].has(pair.snd)) return a;
    const set = a[pair.fst];
    if (set.size === 1) {
        return ofObjectOrdered(OrderedSet.subtract(a.keys, OrderedSet.ofSingleton(pair.fst)), a);
    } else {
        return { ...a, [pair.fst]: OrderedSet.subtract(set, OrderedSet.ofSingleton(pair.snd)), size: a.size - 1, hashCode: -1 }
    }
}

function subtractEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return MultiSet.Empty;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));

    const keys = [], ret = Object.create(null);
    for (let i = 0; i < start; i++) {
        const k = keysA.elementAt(i);
        keys[keys.length] = k;
        ret[k] = a[k];
    }
    for (let i = start; i < end; i++) {
        const k = keysA.elementAt(i);
        if (keysB.has(k)) {
            const subtraction = OrderedSet.subtract(a[k], b[k]);
            if (subtraction.size > 0) {
                keys[keys.length] = k;
                ret[k] = subtraction;
            }
        } else {
            keys[keys.length] = k;
            ret[k] = a[k];
        }
    }
    for (let i = end, _i = keysA.size; i < _i; i++) {
        const k = keysA.elementAt(i);
        keys[keys.length] = k;
        ret[k] = a[k];
    }
    return ofObjectOrdered(OrderedSet.ofSortedArray(keys), ret);
}

export default MultiSet