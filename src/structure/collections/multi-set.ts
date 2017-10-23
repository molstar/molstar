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
        return OrderedSet.has(set.keys, key);
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

    export function union(a: MultiSet, b: MultiSet): MultiSet {
        return findUnion([a, b]);
    }

    export function unionMany(sets: ArrayLike<MultiSet>): MultiSet {
        return findUnion(sets);
    }

    class ElementsIterator implements Iterator<IntPair> {
        private pair = IntPair.zero();
        private unit = 0;

        private keyCount: number;
        private setIndex = -1;
        private currentIndex = 0;
        private currentSize = 0;
        private currentSet: OrderedSet = OrderedSet.Empty;

        [Symbol.iterator]() { return new ElementsIterator(this.elements); };
        done: boolean;
        next() { const value = this.move(); return { value, done: this.done } }

        move() {
            if (this.done) return this.pair;

            if (this.currentIndex >= this.currentSize) {
                if (!this.advance()) return this.pair;
            }

            const next = OrderedSet.elementAt(this.currentSet, this.currentIndex++);
            this.pair.snd = next;
            return this.pair;
        }

        private advance() {
            const keys = this.elements.keys;
            if (++this.setIndex >= this.keyCount) {
                this.done = true;
                return false;
            }
            this.unit = OrderedSet.elementAt(keys, this.setIndex);
            this.pair.fst = this.unit;
            this.currentSet = this.elements[this.unit];
            this.currentIndex = 0;
            this.currentSize = OrderedSet.size(this.currentSet);
            return true;
        }

        constructor(private elements: MultiSetElements) {
            this.keyCount = OrderedSet.size(elements.keys);
            this.done = this.keyCount === 0;
            this.advance();
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
        if (OrderedSet.size(data[k]) > 0) keys[keys.length] = k;
    }
    if (!keys.length) return MultiSet.Empty;
    if (keys.length === 1) {
        const set = data[keys[0]];
        if (OrderedSet.size(set) === 1) return IntPair.pack1(keys[0], OrderedSet.elementAt(set, 0));
    }
    return ofObject1(keys, data);
}

function ofObject1(keys: number[], data: { [id: number]: OrderedSet }) {
    if (keys.length === 1) {
        const k = keys[0];
        const set = data[k];
        if (OrderedSet.size(set) === 1) return IntPair.pack1(k, OrderedSet.elementAt(set, 0));
    }
    sortArray(keys);
    return _createObjectOrdered(OrderedSet.ofSortedArray(keys), data);
}

function ofObjectOrdered(keys: OrderedSet, data: { [id: number]: OrderedSet }) {
    if (OrderedSet.size(keys) === 1) {
        const k = OrderedSet.elementAt(keys, 0);
        const set = data[k];
        if (OrderedSet.size(set) === 1) return IntPair.pack1(k, OrderedSet.elementAt(set, 0));
    }
    return _createObjectOrdered(keys, data);
}

function _createObjectOrdered(keys: OrderedSet, data: { [id: number]: OrderedSet }) {
    const ret: MultiSetElements = Object.create(null);
    ret.keys = keys;
    let size = 0;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.elementAt(keys, i);
        const set = data[k];
        ret[k] = set;
        size += OrderedSet.size(set);
    }
    ret.size = size;
    ret.hashCode = -1;
    return ret;
}

function getUniqueElements(xs: number[]) {
    let count = 1;
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i - 1] !== xs[i]) count++;
    }
    const ret = new (xs as any).constructor(count);
    ret[0] = xs[0];
    let offset = 1;
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i - 1] !== xs[i]) ret[offset++] = xs[i];
    }
    return ret;
}

function normalizeArray(xs: number[]) {
    sortArray(xs);
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i - 1] === xs[i]) return getUniqueElements(xs);
    }
    return xs;
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
        ret[k] = OrderedSet.ofSortedArray(normalizeArray(sets[k]));
    }
    return ofObject1(keys, ret);
}

function computeHash(set: MultiSetElements) {
    const { keys } = set;
    let hash = 23;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.elementAt(keys, i);
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
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.elementAt(keys, i);
        if (!OrderedSet.areEqual(a[k], b[k])) return false;
    }
    return true;
}

function areIntersectingNE(a: number, b: MultiSetElements) {
    IntPair.unpack(a, pair);
    return OrderedSet.has(b.keys, pair.fst) && OrderedSet.has(b[pair.fst], pair.snd);
}

function areIntersectingEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return true;
    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return false;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));
    for (let i = start; i < end; i++) {
        const k = OrderedSet.elementAt(keysA, i);
        if (OrderedSet.has(keysB, k) && OrderedSet.areIntersecting(a[k], b[k])) return true;
    }
    return false;
}

function intersectNE(a: number, b: MultiSetElements) {
    IntPair.unpack(a, pair);
    return OrderedSet.has(b.keys, pair.fst) && OrderedSet.has(b[pair.fst], pair.snd) ? a : MultiSet.Empty;
}

function intersectEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return MultiSet.Empty;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));

    const keys = [], ret = Object.create(null);
    for (let i = start; i < end; i++) {
        const k = OrderedSet.elementAt(keysA, i);
        if (OrderedSet.has(keysB, k)) {
            const intersection = OrderedSet.intersect(a[k], b[k]);
            if (OrderedSet.size(intersection) > 0) {
                keys[keys.length] = k;
                ret[k] = intersection;
            }
        }
    }
    return ofObjectOrdered(OrderedSet.ofSortedArray(keys), ret);
}

function subtractNE(a: number, b: MultiSetElements) {
    IntPair.unpack(a, pair);
    return OrderedSet.has(b.keys, pair.fst) && OrderedSet.has(b[pair.fst], pair.snd) ? MultiSet.Empty : a;
}

function subtractEN(a: MultiSetElements, b: number): MultiSet {
    const aKeys =  a.keys;
    IntPair.unpack(b, pair);
    if (!OrderedSet.has(aKeys, pair.fst) || !OrderedSet.has(a[pair.fst], pair.snd)) return a;
    const set = a[pair.fst];
    if (OrderedSet.size(set) === 1) {
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
        const k = OrderedSet.elementAt(keysA, i);
        keys[keys.length] = k;
        ret[k] = a[k];
    }
    for (let i = start; i < end; i++) {
        const k = OrderedSet.elementAt(keysA, i);
        if (OrderedSet.has(keysB, k)) {
            const subtraction = OrderedSet.subtract(a[k], b[k]);
            if (OrderedSet.size(subtraction) > 0) {
                keys[keys.length] = k;
                ret[k] = subtraction;
            }
        } else {
            keys[keys.length] = k;
            ret[k] = a[k];
        }
    }
    for (let i = end, _i = OrderedSet.size(keysA); i < _i; i++) {
        const k = OrderedSet.elementAt(keysA, i);
        keys[keys.length] = k;
        ret[k] = a[k];
    }
    return ofObjectOrdered(OrderedSet.ofSortedArray(keys), ret);
}

function findUnion(sets: ArrayLike<MultiSet>) {
    if (!sets.length) return MultiSet.Empty;
    if (sets.length === 1) return sets[0];
    if (sets.length === 2 && sets[0] === sets[1]) return sets[0];

    const eCount = { count: 0 };
    const ns = unionN(sets, eCount);
    if (!eCount.count) return ns;
    const ret = Object.create(null);
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s !== 'number') unionInto(ret, s);
    }
    if (MultiSet.size(ns) > 0) {
        if (typeof ns === 'number') unionIntoN(ret, ns);
        else unionInto(ret, ns);
    }
    return ofObject(ret);
}

function unionN(sets: ArrayLike<MultiSet>, eCount: { count: number }) {
    let countN = 0, countE = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        if (typeof sets[i] === 'number') countN++;
        else countE++;
    }
    eCount.count = countE;
    if (!countN) return MultiSet.Empty;
    if (countN === sets.length) return ofPackedPairs(sets as ArrayLike<number>);
    const packed = new Float64Array(countN);
    let offset = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s === 'number') packed[offset++] = s;
    }
    return ofPackedPairs(packed);
}

function unionInto(data: { [key: number]: OrderedSet }, a: MultiSetElements) {
    const keys = a.keys;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.elementAt(keys, i);
        const set = data[k];
        if (set) data[k] = OrderedSet.union(set, a[k]);
        else data[k] = a[k];
    }
}

function unionIntoN(data: { [key: number]: OrderedSet }, a: number) {
    IntPair.unpack(a, pair);
    const set = data[pair.fst];
    if (set) {
        data[pair.fst] = OrderedSet.union(set, OrderedSet.ofSingleton(pair.snd));
    } else {
        data[pair.fst] = OrderedSet.ofSingleton(pair.snd);
    }
}

export default MultiSet