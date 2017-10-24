/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from './ordered-set'
import Iterator from './iterator'
import IntTuple from './int-tuple'
import { sortArray } from './sort'
import { hash1 } from './hash-functions'

/** A map-like representation of integer set */
interface MultiSet { '@type': 'int-multi-set' }

namespace MultiSet {
    export const Empty: MultiSet = { offsets: [0], hashCode: 0, keys: OrderedSet.Empty } as MultiSetElements as any;

    export function create(data: IntTuple | ArrayLike<IntTuple> | IntTuple | { [id: number]: OrderedSet }): MultiSet {
        if (typeof data === 'number') return data;
        if (IntTuple.is(data)) return IntTuple.pack1(data) as any;
        if (isArrayLike(data)) return ofTuples(data) as any;
        return ofObject(data as { [id: number]: OrderedSet }) as any;
    }

    export const keys: (set: MultiSet) => OrderedSet = keysI as any;
    export const keyCount: (set: MultiSet) => number = keyCountI as any;
    export const hasKey: (set: MultiSet, key: number) => boolean = hasKeyI as any;
    export const geyKey: (set: MultiSet, i: number) => number = getKeyI as any;
    export const getByKey: (set: MultiSet, key: number) => OrderedSet = getByKeyI as any;
    export const getByIndex: (set: MultiSet, i: number) => OrderedSet = getByIndexI as any;
    export const has: (set: MultiSet, x: IntTuple) => boolean = hasI as any;
    export const indexOf: (set: MultiSet, x: IntTuple) => number = indexOfI as any;
    export const getAt: (set: MultiSet, i: number) => IntTuple = getAtI as any;
    export const values: (set: MultiSet) => Iterator<IntTuple.Unpacked> = valuesI as any;

    export const size: (set: MultiSet) => number = sizeI as any;
    export const hashCode: (set: MultiSet) => number = hashCodeI as any;

    export const areEqual: (a: MultiSet, b: MultiSet) => boolean = areEqualI as any;
    export const areIntersecting: (a: MultiSet, b: MultiSet) => boolean = areIntersectingI as any;

    export const union: (a: MultiSet, b: MultiSet) => MultiSet = unionI as any;
    export const unionMany: (sets: MultiSet[]) => MultiSet = unionManyI as any;
    export const intersect: (a: MultiSet, b: MultiSet) => MultiSet = intersectI as any;
    export const subtract: (a: MultiSet, b: MultiSet) => MultiSet = subtractI as any;
}

export default MultiSet

/** Long and painful implementation starts here */

interface MultiSetElements { [id: number]: OrderedSet, offsets: number[], hashCode: number, keys: OrderedSet }
type MultiSetImpl = IntTuple | MultiSetElements

function keysI(set: MultiSetImpl): OrderedSet {
    if (typeof set === 'number') return OrderedSet.ofSingleton(set);
    return (set as MultiSetElements).keys;
}

function keyCountI(set: MultiSetImpl): number {
    if (typeof set === 'number') return 1;
    return OrderedSet.size((set as MultiSetElements).keys);
}

function hasKeyI(set: MultiSetImpl, key: number): boolean {
    if (typeof set === 'number') return IntTuple.fst(set) === key;
    return OrderedSet.has((set as MultiSetElements).keys, key);
}

function getKeyI(set: MultiSetImpl, index: number): number {
    if (typeof set === 'number') return IntTuple.fst(set);
    return OrderedSet.getAt((set as MultiSetElements).keys, index);
}

function hasI(set: MultiSetImpl, t: IntTuple): boolean {
    if (typeof set === 'number') return IntTuple.areEqual(t, set);
    IntTuple.unpack(t, _hasP);
    return OrderedSet.has((set as MultiSetElements).keys, _hasP.fst) ? OrderedSet.has((set as MultiSetElements)[_hasP.fst], _hasP.snd) : false;
}
const _hasP = IntTuple.zero();

function getByKeyI(set: MultiSetImpl, key: number): OrderedSet {
    if (typeof set === 'number') {
        IntTuple.unpack(set, _gS);
        return _gS.fst === key ? OrderedSet.ofSingleton(_gS.snd) : OrderedSet.Empty;
    }
    return OrderedSet.has((set as MultiSetElements).keys, key) ? (set as MultiSetElements)[key] : OrderedSet.Empty;
}
const _gS = IntTuple.zero();

function getByIndexI(set: MultiSetImpl, index: number): OrderedSet {
    if (typeof set === 'number') return index === 0 ? OrderedSet.ofSingleton(IntTuple.snd(set)) : OrderedSet.Empty;
    const key = OrderedSet.getAt((set as MultiSetElements).keys, index);
    return (set as MultiSetElements)[key] || OrderedSet.Empty;
}

function getAtI(set: MultiSetImpl, i: number): IntTuple {
    if (typeof set === 'number') return set;
    return getAtE(set as MultiSetElements, i);
}

function indexOfI(set: MultiSetImpl, t: IntTuple) {
    if (typeof set === 'number') return IntTuple.areEqual(set, t) ? 0 : -1;
    return indexOfE(set as MultiSetElements, t);
}

/** Number elements in the "child" sets */
function sizeI(set: MultiSetImpl) {
    if (typeof set === 'number') return 1;
    return (set as MultiSetElements).offsets[(set as MultiSetElements).offsets.length - 1];
}

function hashCodeI(set: MultiSetImpl) {
    if (typeof set === 'number') return IntTuple.hashCode(set);
    if ((set as MultiSetElements).hashCode !== -1) return (set as MultiSetElements).hashCode;
    return computeHash((set as MultiSetElements));
}

function areEqualI(a: MultiSetImpl, b: MultiSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b);
        return false;
    }
    if (typeof b === 'number') return false;
    return areEqualEE(a as MultiSetElements, b as MultiSetElements);
}

function areIntersectingI(a: MultiSetImpl, b: MultiSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b);
        return areIntersectingNE(a, b as MultiSetElements);
    }
    if (typeof b === 'number') return areIntersectingNE(b, a as MultiSetElements);
    return areIntersectingEE(a as MultiSetElements, b as MultiSetElements);
}

function intersectI(a: MultiSetImpl, b: MultiSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b) ? a : MultiSet.Empty;
        return intersectNE(a, b as MultiSetElements);
    }
    if (typeof b === 'number') return intersectNE(b, a as MultiSetElements);
    return intersectEE(a as MultiSetElements, b as MultiSetElements);
}

function subtractI(a: MultiSetImpl, b: MultiSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b) ? MultiSet.Empty : a;
        return subtractNE(a, b as MultiSetElements);
    }
    if (typeof b === 'number') return subtractEN(a as MultiSetElements, b);
    return subtractEE(a as MultiSetElements, b as MultiSetElements);
}

function unionI(a: MultiSetImpl, b: MultiSetImpl) {
    return findUnion([a, b]);
}

function unionManyI(sets: ArrayLike<MultiSetImpl>) {
    return findUnion(sets);
}

class ElementsIterator implements Iterator<IntTuple.Unpacked> {
    private pair = IntTuple.zero();

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

        this.pair.snd = OrderedSet.getAt(this.currentSet, this.currentIndex++);
        return this.pair;
    }

    private advance() {
        if (++this.setIndex >= this.keyCount) {
            this.done = true;
            return false;
        }
        const unit = OrderedSet.getAt(this.elements.keys, this.setIndex);
        this.pair.fst = unit;
        this.currentSet = this.elements[unit];
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

function valuesI(set: MultiSetImpl): Iterator<IntTuple.Unpacked> {
    if (typeof set === 'number') return Iterator.Value(IntTuple.unpack1(set));
    return new ElementsIterator(set as MultiSetElements);
}

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
        if (OrderedSet.size(set) === 1) return IntTuple.pack(keys[0], OrderedSet.getAt(set, 0));
    }
    return ofObject1(keys, data);
}

function ofObject1(keys: number[], data: { [id: number]: OrderedSet }) {
    if (keys.length === 1) {
        const k = keys[0];
        const set = data[k];
        if (OrderedSet.size(set) === 1) return IntTuple.pack(k, OrderedSet.getAt(set, 0));
    }
    sortArray(keys);
    return _createObjectOrdered(OrderedSet.ofSortedArray(keys), data);
}

function ofObjectOrdered(keys: OrderedSet, data: { [id: number]: OrderedSet }) {
    if (OrderedSet.size(keys) === 1) {
        const k = OrderedSet.getAt(keys, 0);
        const set = data[k];
        if (OrderedSet.size(set) === 1) return IntTuple.pack(k, OrderedSet.getAt(set, 0));
    }
    return _createObjectOrdered(keys, data);
}

function _createObjectOrdered(keys: OrderedSet, data: { [id: number]: OrderedSet }) {
    const ret: MultiSetElements = Object.create(null);
    ret.keys = keys;
    const offsets = [0];
    let size = 0;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.getAt(keys, i);
        const set = data[k];
        ret[k] = set;
        size += OrderedSet.size(set);
        offsets[offsets.length] = size;
    }
    ret.offsets = offsets;
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

function ofTuples(xs: ArrayLike<IntTuple>) {
    if (xs.length === 0) return MultiSet.Empty;
    const sets: { [key: number]: number[] } = Object.create(null);
    const p = IntTuple.zero();
    for (let i = 0, _i = xs.length; i < _i; i++) {
        IntTuple.unpack(xs[i], p);
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

function getOffsetIndex(xs: ArrayLike<number>, value: number) {
    let min = 0, max = xs.length - 1;
    while (min < max) {
        const mid = (min + max) >> 1;
        const v = xs[mid];
        if (value < v) max = mid - 1;
        else if (value > v) min = mid + 1;
        else return mid;
    }
    if (min > max) {
        return max;
    }
    return value < xs[min] ? min - 1 : min;
}

function getAtE(set: MultiSetElements, i: number) {
    const { offsets, keys } = set;
    const o = getOffsetIndex(offsets, i);
    if (o >= offsets.length - 1) return 0;
    const k = OrderedSet.getAt(keys, o);
    const e = OrderedSet.getAt(set[k], i - offsets[o]);
    return IntTuple.pack(k, e);
}

const _iOE = IntTuple.zero();
function indexOfE(set: MultiSetElements, t: IntTuple) {
    IntTuple.unpack(t, _iOE);
    const { keys } = set;
    const setIdx = OrderedSet.indexOf(keys, _iOE.fst);
    if (setIdx < 0) return -1;
    const o = OrderedSet.indexOf(set[_iOE.fst], _iOE.snd);
    if (o < 0) return -1;
    return set.offsets[setIdx] + o;
}

function computeHash(set: MultiSetElements) {
    const { keys } = set;
    let hash = 23;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.getAt(keys, i);
        hash = (31 * hash + k) | 0;
        hash = (31 * hash + OrderedSet.hashCode(set[k])) | 0;
    }
    hash = (31 * hash + sizeI(set)) | 0;
    hash = hash1(hash);
    set.hashCode = hash;
    return hash;
}

function areEqualEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return true;
    if (sizeI(a) !== sizeI(a)) return false;

    const keys = a.keys;
    if (!OrderedSet.areEqual(keys, b.keys)) return false;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.getAt(keys, i);
        if (!OrderedSet.areEqual(a[k], b[k])) return false;
    }
    return true;
}

const _aeP = IntTuple.zero();
function areIntersectingNE(a: IntTuple, b: MultiSetElements) {
    IntTuple.unpack(a, _aeP);
    return OrderedSet.has(b.keys, _aeP.fst) && OrderedSet.has(b[_aeP.fst], _aeP.snd);
}

function areIntersectingEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return true;
    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return false;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));
    for (let i = start; i < end; i++) {
        const k = OrderedSet.getAt(keysA, i);
        if (OrderedSet.has(keysB, k) && OrderedSet.areIntersecting(a[k], b[k])) return true;
    }
    return false;
}

const _nP = IntTuple.zero();
function intersectNE(a: IntTuple, b: MultiSetElements) {
    IntTuple.unpack(a, _nP);
    return OrderedSet.has(b.keys, _nP.fst) && OrderedSet.has(b[_nP.fst], _nP.snd) ? a : MultiSet.Empty;
}

function intersectEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return MultiSet.Empty;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));

    const keys = [], ret = Object.create(null);
    for (let i = start; i < end; i++) {
        const k = OrderedSet.getAt(keysA, i);
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

const _sNE = IntTuple.zero();
function subtractNE(a: IntTuple, b: MultiSetElements) {
    IntTuple.unpack(a, _sNE);
    return OrderedSet.has(b.keys, _sNE.fst) && OrderedSet.has(b[_sNE.fst], _sNE.snd) ? MultiSet.Empty : a;
}

const _sEN = IntTuple.zero();
function subtractEN(a: MultiSetElements, b: IntTuple): MultiSetImpl {
    const aKeys =  a.keys;
    IntTuple.unpack(b, _sEN);
    if (!OrderedSet.has(aKeys, _sEN.fst) || !OrderedSet.has(a[_sEN.fst], _sEN.snd)) return a;
    const set = a[_sEN.fst];
    if (OrderedSet.size(set) === 1) {
        return ofObjectOrdered(OrderedSet.subtract(a.keys, OrderedSet.ofSingleton(_sEN.fst)), a);
    } else {
        return ofObjectOrdered(OrderedSet.subtract(set, OrderedSet.ofSingleton(_sEN.snd)), a);
    }
}

function subtractEE(a: MultiSetElements, b: MultiSetElements) {
    if (a === b) return MultiSet.Empty;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return MultiSet.Empty;
    const { start, end } = OrderedSet.getIntervalRange(keysA, OrderedSet.min(keysB), OrderedSet.max(keysB));

    const keys = [], ret = Object.create(null);
    for (let i = 0; i < start; i++) {
        const k = OrderedSet.getAt(keysA, i);
        keys[keys.length] = k;
        ret[k] = a[k];
    }
    for (let i = start; i < end; i++) {
        const k = OrderedSet.getAt(keysA, i);
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
        const k = OrderedSet.getAt(keysA, i);
        keys[keys.length] = k;
        ret[k] = a[k];
    }
    return ofObjectOrdered(OrderedSet.ofSortedArray(keys), ret);
}

function findUnion(sets: ArrayLike<MultiSetImpl>) {
    if (!sets.length) return MultiSet.Empty;
    if (sets.length === 1) return sets[0];
    if (sets.length === 2 && areEqualI(sets[0], sets[1])) return sets[0];

    const eCount = { count: 0 };
    const ns = unionN(sets, eCount);
    if (!eCount.count) return ns;
    const ret = Object.create(null);
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s !== 'number') unionInto(ret, s as MultiSetElements);
    }
    if (sizeI(ns as MultiSetImpl) > 0) {
        if (typeof ns === 'number') unionIntoN(ret, ns);
        else unionInto(ret, ns as MultiSetElements);
    }
    return ofObject(ret);
}

function unionN(sets: ArrayLike<MultiSetImpl>, eCount: { count: number }) {
    let countN = 0, countE = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        if (typeof sets[i] === 'number') countN++;
        else countE++;
    }
    eCount.count = countE;
    if (!countN) return MultiSet.Empty;
    if (countN === sets.length) return ofTuples(sets as ArrayLike<number>);
    const packed = new Float64Array(countN);
    let offset = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s === 'number') packed[offset++] = s;
    }
    return ofTuples(packed);
}

function unionInto(data: { [key: number]: OrderedSet }, a: MultiSetElements) {
    const keys = a.keys;
    for (let i = 0, _i = OrderedSet.size(keys); i < _i; i++) {
        const k = OrderedSet.getAt(keys, i);
        const set = data[k];
        if (set) data[k] = OrderedSet.union(set, a[k]);
        else data[k] = a[k];
    }
}

const _uIN = IntTuple.zero();
function unionIntoN(data: { [key: number]: OrderedSet }, a: IntTuple) {
    IntTuple.unpack(a, _uIN);
    const set = data[_uIN.fst];
    if (set) {
        data[_uIN.fst] = OrderedSet.union(set, OrderedSet.ofSingleton(_uIN.snd));
    } else {
        data[_uIN.fst] = OrderedSet.ofSingleton(_uIN.snd);
    }
}