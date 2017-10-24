/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from '../common/collections/ordered-set'
import Iterator from '../common/collections/iterator'
import IntTuple from '../common/collections/int-tuple'
import { sortArray } from '../common/collections/sort'
import { hash1 } from '../common/collections/hash-functions'

/** A map-like representation of integer set */
interface AtomSet { '@type': 'atom-set' }

namespace AtomSet {
    export const Empty: AtomSet = { offsets: [0], hashCode: 0, keys: OrderedSet.Empty } as AtomSetElements as any;

    export function create(data: IntTuple | ArrayLike<IntTuple> | IntTuple | { [id: number]: OrderedSet }): AtomSet {
        if (typeof data === 'number') return data;
        if (IntTuple.is(data)) return IntTuple.pack1(data) as any;
        if (isArrayLike(data)) return ofTuples(data) as any;
        return ofObject(data as { [id: number]: OrderedSet }) as any;
    }

    export const keys: (set: AtomSet) => OrderedSet = keysI as any;
    export const keyCount: (set: AtomSet) => number = keyCountI as any;
    export const hasKey: (set: AtomSet, key: number) => boolean = hasKeyI as any;
    export const geyKey: (set: AtomSet, i: number) => number = getKeyI as any;
    export const getByKey: (set: AtomSet, key: number) => OrderedSet = getByKeyI as any;
    export const getByIndex: (set: AtomSet, i: number) => OrderedSet = getByIndexI as any;
    export const has: (set: AtomSet, x: IntTuple) => boolean = hasI as any;
    export const indexOf: (set: AtomSet, x: IntTuple) => number = indexOfI as any;
    export const getAt: (set: AtomSet, i: number) => IntTuple = getAtI as any;
    export const values: (set: AtomSet) => Iterator<IntTuple.Unpacked> = valuesI as any;

    export const size: (set: AtomSet) => number = sizeI as any;
    export const hashCode: (set: AtomSet) => number = hashCodeI as any;

    export const areEqual: (a: AtomSet, b: AtomSet) => boolean = areEqualI as any;
    export const areIntersecting: (a: AtomSet, b: AtomSet) => boolean = areIntersectingI as any;

    export const union: (a: AtomSet, b: AtomSet) => AtomSet = unionI as any;
    export const unionMany: (sets: AtomSet[]) => AtomSet = unionManyI as any;
    export const intersect: (a: AtomSet, b: AtomSet) => AtomSet = intersectI as any;
    export const subtract: (a: AtomSet, b: AtomSet) => AtomSet = subtractI as any;
}

export default AtomSet

/** Long and painful implementation starts here */

interface AtomSetElements { [id: number]: OrderedSet, offsets: number[], hashCode: number, keys: OrderedSet }
type AtomSetImpl = IntTuple | AtomSetElements

function keysI(set: AtomSetImpl): OrderedSet {
    if (typeof set === 'number') return OrderedSet.ofSingleton(set);
    return (set as AtomSetElements).keys;
}

function keyCountI(set: AtomSetImpl): number {
    if (typeof set === 'number') return 1;
    return OrderedSet.size((set as AtomSetElements).keys);
}

function hasKeyI(set: AtomSetImpl, key: number): boolean {
    if (typeof set === 'number') return IntTuple.fst(set) === key;
    return OrderedSet.has((set as AtomSetElements).keys, key);
}

function getKeyI(set: AtomSetImpl, index: number): number {
    if (typeof set === 'number') return IntTuple.fst(set);
    return OrderedSet.getAt((set as AtomSetElements).keys, index);
}

function hasI(set: AtomSetImpl, t: IntTuple): boolean {
    if (typeof set === 'number') return IntTuple.areEqual(t, set);
    IntTuple.unpack(t, _hasP);
    return OrderedSet.has((set as AtomSetElements).keys, _hasP.fst) ? OrderedSet.has((set as AtomSetElements)[_hasP.fst], _hasP.snd) : false;
}
const _hasP = IntTuple.zero();

function getByKeyI(set: AtomSetImpl, key: number): OrderedSet {
    if (typeof set === 'number') {
        IntTuple.unpack(set, _gS);
        return _gS.fst === key ? OrderedSet.ofSingleton(_gS.snd) : OrderedSet.Empty;
    }
    return OrderedSet.has((set as AtomSetElements).keys, key) ? (set as AtomSetElements)[key] : OrderedSet.Empty;
}
const _gS = IntTuple.zero();

function getByIndexI(set: AtomSetImpl, index: number): OrderedSet {
    if (typeof set === 'number') return index === 0 ? OrderedSet.ofSingleton(IntTuple.snd(set)) : OrderedSet.Empty;
    const key = OrderedSet.getAt((set as AtomSetElements).keys, index);
    return (set as AtomSetElements)[key] || OrderedSet.Empty;
}

function getAtI(set: AtomSetImpl, i: number): IntTuple {
    if (typeof set === 'number') return set;
    return getAtE(set as AtomSetElements, i);
}

function indexOfI(set: AtomSetImpl, t: IntTuple) {
    if (typeof set === 'number') return IntTuple.areEqual(set, t) ? 0 : -1;
    return indexOfE(set as AtomSetElements, t);
}

/** Number elements in the "child" sets */
function sizeI(set: AtomSetImpl) {
    if (typeof set === 'number') return 1;
    return (set as AtomSetElements).offsets[(set as AtomSetElements).offsets.length - 1];
}

function hashCodeI(set: AtomSetImpl) {
    if (typeof set === 'number') return IntTuple.hashCode(set);
    if ((set as AtomSetElements).hashCode !== -1) return (set as AtomSetElements).hashCode;
    return computeHash((set as AtomSetElements));
}

function areEqualI(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b);
        return false;
    }
    if (typeof b === 'number') return false;
    return areEqualEE(a as AtomSetElements, b as AtomSetElements);
}

function areIntersectingI(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b);
        return areIntersectingNE(a, b as AtomSetElements);
    }
    if (typeof b === 'number') return areIntersectingNE(b, a as AtomSetElements);
    return areIntersectingEE(a as AtomSetElements, b as AtomSetElements);
}

function intersectI(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b) ? a : AtomSet.Empty;
        return intersectNE(a, b as AtomSetElements);
    }
    if (typeof b === 'number') return intersectNE(b, a as AtomSetElements);
    return intersectEE(a as AtomSetElements, b as AtomSetElements);
}

function subtractI(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return IntTuple.areEqual(a, b) ? AtomSet.Empty : a;
        return subtractNE(a, b as AtomSetElements);
    }
    if (typeof b === 'number') return subtractEN(a as AtomSetElements, b);
    return subtractEE(a as AtomSetElements, b as AtomSetElements);
}

function unionI(a: AtomSetImpl, b: AtomSetImpl) {
    return findUnion([a, b]);
}

function unionManyI(sets: ArrayLike<AtomSetImpl>) {
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

    constructor(private elements: AtomSetElements) {
        this.keyCount = OrderedSet.size(elements.keys);
        this.done = this.keyCount === 0;
        this.advance();
    }
}

function valuesI(set: AtomSetImpl): Iterator<IntTuple.Unpacked> {
    if (typeof set === 'number') return Iterator.Value(IntTuple.unpack1(set));
    return new ElementsIterator(set as AtomSetElements);
}

function isArrayLike(x: any): x is ArrayLike<IntTuple> {
    return x && (typeof x.length === 'number' && (x instanceof Array || !!x.buffer));
}

function ofObject(data: { [id: number]: OrderedSet }) {
    const keys = [];
    for (const _k of Object.keys(data)) {
        const k = +_k;
        if (OrderedSet.size(data[k]) > 0) keys[keys.length] = k;
    }
    if (!keys.length) return AtomSet.Empty;
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
    const ret: AtomSetElements = Object.create(null);
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
    if (xs.length === 0) return AtomSet.Empty;
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

function getAtE(set: AtomSetElements, i: number): IntTuple {
    const { offsets, keys } = set;
    const o = getOffsetIndex(offsets, i);
    if (o >= offsets.length - 1) return 0 as any;
    const k = OrderedSet.getAt(keys, o);
    const e = OrderedSet.getAt(set[k], i - offsets[o]);
    return IntTuple.pack(k, e);
}

const _iOE = IntTuple.zero();
function indexOfE(set: AtomSetElements, t: IntTuple) {
    IntTuple.unpack(t, _iOE);
    const { keys } = set;
    const setIdx = OrderedSet.indexOf(keys, _iOE.fst);
    if (setIdx < 0) return -1;
    const o = OrderedSet.indexOf(set[_iOE.fst], _iOE.snd);
    if (o < 0) return -1;
    return set.offsets[setIdx] + o;
}

function computeHash(set: AtomSetElements) {
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

function areEqualEE(a: AtomSetElements, b: AtomSetElements) {
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
function areIntersectingNE(a: IntTuple, b: AtomSetElements) {
    IntTuple.unpack(a, _aeP);
    return OrderedSet.has(b.keys, _aeP.fst) && OrderedSet.has(b[_aeP.fst], _aeP.snd);
}

function areIntersectingEE(a: AtomSetElements, b: AtomSetElements) {
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
function intersectNE(a: IntTuple, b: AtomSetElements) {
    IntTuple.unpack(a, _nP);
    return OrderedSet.has(b.keys, _nP.fst) && OrderedSet.has(b[_nP.fst], _nP.snd) ? a : AtomSet.Empty;
}

function intersectEE(a: AtomSetElements, b: AtomSetElements) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return AtomSet.Empty;
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
function subtractNE(a: IntTuple, b: AtomSetElements) {
    IntTuple.unpack(a, _sNE);
    return OrderedSet.has(b.keys, _sNE.fst) && OrderedSet.has(b[_sNE.fst], _sNE.snd) ? AtomSet.Empty : a;
}

const _sEN = IntTuple.zero();
function subtractEN(a: AtomSetElements, b: IntTuple): AtomSetImpl {
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

function subtractEE(a: AtomSetElements, b: AtomSetElements) {
    if (a === b) return AtomSet.Empty;

    const keysA = a.keys, keysB = b.keys;
    if (!OrderedSet.areIntersecting(a.keys, b.keys)) return AtomSet.Empty;
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

function findUnion(sets: ArrayLike<AtomSetImpl>) {
    if (!sets.length) return AtomSet.Empty;
    if (sets.length === 1) return sets[0];
    if (sets.length === 2 && areEqualI(sets[0], sets[1])) return sets[0];

    const eCount = { count: 0 };
    const ns = unionN(sets, eCount);
    if (!eCount.count) return ns;
    const ret = Object.create(null);
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s !== 'number') unionInto(ret, s as AtomSetElements);
    }
    if (sizeI(ns as AtomSetImpl) > 0) {
        if (typeof ns === 'number') unionIntoN(ret, ns as any);
        else unionInto(ret, ns as AtomSetElements);
    }
    return ofObject(ret);
}

function unionN(sets: ArrayLike<AtomSetImpl>, eCount: { count: number }) {
    let countN = 0, countE = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        if (typeof sets[i] === 'number') countN++;
        else countE++;
    }
    eCount.count = countE;
    if (!countN) return AtomSet.Empty;
    if (countN === sets.length) return ofTuples(sets as ArrayLike<IntTuple>);
    const packed = new Float64Array(countN);
    let offset = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s === 'number') packed[offset++] = s;
    }
    return ofTuples(packed as any);
}

function unionInto(data: { [key: number]: OrderedSet }, a: AtomSetElements) {
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