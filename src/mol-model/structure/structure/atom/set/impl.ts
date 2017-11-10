/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray, Interval, Iterator, OrderedSet, IntMap } from 'mol-data/int'
import { sortArray } from 'mol-data/util/sort'
import { hash1 } from 'mol-data/util/hash-functions'
import Atom from '../../atom'

/** Long and painful implementation starts here */

export interface AtomSetElements { sets: IntMap<OrderedSet>, offsets: Int32Array, hashCode: number, keys: SortedArray }
export type AtomSetImpl = Atom | AtomSetElements

export const Empty: AtomSetImpl = { sets: IntMap.Empty, offsets: new Int32Array(1), hashCode: 0, keys: SortedArray.Empty };

export function create(data: Atom | ArrayLike<Atom>): AtomSetImpl {
    if (typeof data === 'number' || Atom.is(data)) return data;
    return ofAtoms(data);
}

export function isSingleton(set: AtomSetImpl) {
    return typeof set === 'number';
}

export function getKeys(set: AtomSetImpl): SortedArray {
    if (typeof set === 'number') return SortedArray.ofSingleton(set);
    return (set as AtomSetElements).keys;
}

export function keyCount(set: AtomSetImpl): number {
    if (typeof set === 'number') return 1;
    return (set as AtomSetElements).keys.length;
}

export function hasKey(set: AtomSetImpl, key: number): boolean {
    if (typeof set === 'number') return Atom.unit(set) === key;
    return !!(set as AtomSetElements).sets.has(key);
}

export function getKey(set: AtomSetImpl, index: number): number {
    if (typeof set === 'number') return Atom.unit(set);
    return (set as AtomSetElements).keys[index];
}

export function hasAtom(set: AtomSetImpl, t: Atom): boolean {
    if (typeof set === 'number') return Atom.areEqual(t, set);
    const os = (set as AtomSetElements).sets.get(Atom.unit(t));
    return !!os && OrderedSet.has(os, Atom.index(t));
}

export function getByKey(set: AtomSetImpl, key: number): OrderedSet {
    if (typeof set === 'number') {
        return Atom.unit(set) === key ? OrderedSet.ofSingleton(Atom.index(set)) : OrderedSet.Empty;
    }
    return (set as AtomSetElements).sets.get(key) || OrderedSet.Empty;
}

export function getByIndex(set: AtomSetImpl, index: number): OrderedSet {
    if (typeof set === 'number') return index === 0 ? OrderedSet.ofSingleton(Atom.index(set)) : OrderedSet.Empty;
    const key = (set as AtomSetElements).keys[index];
    return (set as AtomSetElements).sets.get(key) || OrderedSet.Empty;
}

export function getAt(set: AtomSetImpl, i: number): Atom {
    if (typeof set === 'number') return set;
    return getAtE(set as AtomSetElements, i);
}

export function indexOf(set: AtomSetImpl, t: Atom) {
    if (typeof set === 'number') return Atom.areEqual(set, t) ? 0 : -1;
    return indexOfE(set as AtomSetElements, t);
}

/** Number elements in the "child" sets */
export function size(set: AtomSetImpl) {
    if (typeof set === 'number') return 1;
    return (set as AtomSetElements).offsets[(set as AtomSetElements).offsets.length - 1];
}

export function hashCode(set: AtomSetImpl) {
    if (typeof set === 'number') return Atom.hashCode(set);
    if ((set as AtomSetElements).hashCode !== -1) return (set as AtomSetElements).hashCode;
    return computeHash((set as AtomSetElements));
}

export function areEqual(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return Atom.areEqual(a, b);
        return false;
    }
    if (typeof b === 'number') return false;
    return areEqualEE(a as AtomSetElements, b as AtomSetElements);
}

export function areIntersecting(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return Atom.areEqual(a, b);
        return areIntersectingNE(a, b as AtomSetElements);
    }
    if (typeof b === 'number') return areIntersectingNE(b, a as AtomSetElements);
    return areIntersectingEE(a as AtomSetElements, b as AtomSetElements);
}

export function intersect(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return Atom.areEqual(a, b) ? a : Empty;
        return intersectNE(a, b as AtomSetElements);
    }
    if (typeof b === 'number') return intersectNE(b, a as AtomSetElements);
    return intersectEE(a as AtomSetElements, b as AtomSetElements);
}

export function subtract(a: AtomSetImpl, b: AtomSetImpl) {
    if (typeof a === 'number') {
        if (typeof b === 'number') return Atom.areEqual(a, b) ? Empty : a;
        return subtractNE(a, b as AtomSetElements);
    }
    if (typeof b === 'number') return subtractEN(a as AtomSetElements, b);
    return subtractEE(a as AtomSetElements, b as AtomSetElements);
}

export function union(a: AtomSetImpl, b: AtomSetImpl) {
    return findUnion([a, b]);
}

export function unionMany(sets: ArrayLike<AtomSetImpl>) {
    return findUnion(sets);
}

class ElementsIterator implements Iterator<Atom> {
    private unit: number = 0;
    private keyCount: number;
    private setIndex = -1;
    private currentIndex = 0;
    private currentSize = 0;
    private currentSet: OrderedSet = OrderedSet.Empty;

    hasNext: boolean = false;

    move() {
        if (!this.hasNext) return Atom.Zero;
        const ret = Atom.create(this.unit, OrderedSet.getAt(this.currentSet, this.currentIndex++));
        if (this.currentIndex >= this.currentSize) this.advance();
        return ret;
    }

    private advance() {
        if (++this.setIndex >= this.keyCount) {
            this.hasNext = false;
            return false;
        }
        this.unit = this.elements.keys[this.setIndex];
        this.currentSet = this.elements.sets.get(this.unit);
        this.currentIndex = 0;
        this.currentSize = OrderedSet.size(this.currentSet);
        return true;
    }

    constructor(private elements: AtomSetElements) {
        this.keyCount = elements.keys.length;
        this.hasNext = this.keyCount > 0;
        this.advance();
    }
}

export function values(set: AtomSetImpl): Iterator<Atom> {
    if (typeof set === 'number') return Iterator.Value(set as Atom);
    return new ElementsIterator(set as AtomSetElements);
}

export class AtomSetGenerator {
    private keys: number[] = [];
    private sets = IntMap.Mutable<OrderedSet>();

    add(unit: number, set: OrderedSet) {
        if (OrderedSet.size(set) === 0) return;
        this.keys[this.keys.length] = unit;
        this.sets.set(unit, set);
    }

    addUnion(unit: number, set: OrderedSet) {
        if (OrderedSet.size(set) === 0) return;

        if (this.sets.has(unit)) {
            this.sets.set(unit, OrderedSet.union(this.sets.get(unit), set));
        } else {
            this.keys[this.keys.length] = unit;
            this.sets.set(unit, set);
        }
    }

    getSet(): AtomSetImpl {
        return ofKeysAndSets(this.keys, this.sets);
    }
}

export function Generator() {
    return new AtomSetGenerator();
}

function ofKeysAndSetsElemements(keys: number[], sets: IntMap<OrderedSet>): AtomSetElements {
    sortArray(keys);
    let runningSize = 0;
    const offsets = new Int32Array(keys.length + 1);
    for (let i = 0, _i = keys.length; i < _i; i++) {
        runningSize += OrderedSet.size(sets.get(keys[i]));
        offsets[i + 1] = runningSize;
    }
    return { keys: SortedArray.ofSortedArray(keys), sets: IntMap.asImmutable(sets), offsets, hashCode: -1 };
}


function ofKeysAndSets(keys: number[], sets: IntMap<OrderedSet>) {
    if (keys.length === 1) {
        const set = sets.get(keys[0]);
        if (OrderedSet.size(set) === 1) return Atom.create(keys[0], OrderedSet.getAt(set, 0));
    }
    return ofKeysAndSetsElemements(keys, sets);
}

function getUniqueElements(xs: number[]): number[] {
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

function normalizeArray(xs: number[]): number[] {
    sortArray(xs);
    for (let i = 1, _i = xs.length; i < _i; i++) {
        if (xs[i - 1] === xs[i]) return getUniqueElements(xs);
    }
    return xs;
}

function ofAtoms(xs: ArrayLike<Atom>) {
    if (xs.length === 0) return Empty;

    const elements = IntMap.Mutable<number[]>();
    const keys: number[] = [];
    for (let i = 0, _i = xs.length; i < _i; i++) {
        const x = xs[i];
        const u = Atom.unit(x), v = Atom.index(x);
        if (elements.has(u)) {
            const set = elements.get(u);
            set[set.length] = v;
        } else {
            keys[keys.length] = u;
            elements.set(u, [v]);
        }
    }

    const sets = IntMap.Mutable<OrderedSet>();
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        sets.set(k, OrderedSet.ofSortedArray(normalizeArray(elements.get(k))));
    }

    return ofKeysAndSets(keys, sets);
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

function getAtE(set: AtomSetElements, i: number): Atom {
    const { offsets, keys } = set;
    const o = getOffsetIndex(offsets, i);
    if (o >= offsets.length - 1) return 0 as any;
    const k = keys[o];
    const e = OrderedSet.getAt(set.sets.get(k), i - offsets[o]);
    return Atom.create(k, e);
}

function indexOfE(set: AtomSetElements, t: Atom) {
    const { keys } = set;
    const u = Atom.unit(t);
    const setIdx = SortedArray.indexOf(keys, u);
    if (setIdx < 0) return -1;
    const o = OrderedSet.indexOf(set.sets.get(u), Atom.index(t));
    if (o < 0) return -1;
    return set.offsets[setIdx] + o;
}

function computeHash(set: AtomSetElements) {
    const { keys, sets } = set;
    let hash = 23;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        hash = (31 * hash + k) | 0;
        hash = (31 * hash + OrderedSet.hashCode(sets.get(k))) | 0;
    }
    hash = (31 * hash + size(set)) | 0;
    hash = hash1(hash);
    set.hashCode = hash;
    return hash;
}

function areEqualEE(a: AtomSetElements, b: AtomSetElements) {
    if (a === b) return true;
    if (size(a) !== size(a)) return false;

    const keys = a.keys;
    if (!SortedArray.areEqual(keys, b.keys)) return false;
    const { sets: aSets } = a;
    const { sets: bSets } = b;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        if (!OrderedSet.areEqual(aSets.get(k), bSets.get(k))) return false;
    }
    return true;
}

function areIntersectingNE(a: Atom, b: AtomSetElements) {
    const u = Atom.unit(a);
    return b.sets.has(u) && OrderedSet.has(b.sets.get(u), Atom.index(a));
}

function areIntersectingEE(a: AtomSetElements, b: AtomSetElements) {
    if (a === b) return true;
    const keysA = a.keys, keysB = b.keys;
    if (!SortedArray.areIntersecting(a.keys, b.keys)) return false;
    const r = SortedArray.findRange(keysA, SortedArray.min(keysB), SortedArray.max(keysB));
    const start = Interval.start(r), end = Interval.end(r);
    const { sets: aSets } = a;
    const { sets: bSets } = b;
    for (let i = start; i < end; i++) {
        const k = keysA[i];
        const ak = aSets.get(k), bk = bSets.get(k);
        if (!!ak && !!bk && OrderedSet.areIntersecting(ak, bk)) return true;
    }
    return false;
}

function intersectNE(a: Atom, b: AtomSetElements) {
    const u = Atom.unit(a);
    return b.sets.has(u) && OrderedSet.has(b.sets.get(u), Atom.index(a)) ? a : Empty;
}

function intersectEE(a: AtomSetElements, b: AtomSetElements) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!SortedArray.areIntersecting(a.keys, b.keys)) return Empty;
    const r = SortedArray.findRange(keysA, SortedArray.min(keysB), SortedArray.max(keysB));
    const start = Interval.start(r), end = Interval.end(r);

    const { sets: aSets } = a;
    const { sets: bSets } = b;
    const generator = Generator();
    for (let i = start; i < end; i++) {
        const k = keysA[i];
        const bk = bSets.get(k);
        if (!bk) continue;
        generator.add(k, OrderedSet.intersect(aSets.get(k), bk));
    }
    return generator.getSet();
}

function subtractNE(a: Atom, b: AtomSetElements) {
    return hasAtom(b, a) ? Empty : a;
}

function subtractEN(a: AtomSetElements, b: Atom): AtomSetImpl {
    if (!hasAtom(a, b)) return a;

    const u = Atom.unit(b), v = Atom.index(b);
    const { sets: aSets } = a;
    const set = aSets.get(u);

    if (OrderedSet.size(set) === 1) {
        // remove the entire unit.
        const generator = Generator();
        for (let i = 0, _i = a.keys.length; i < _i; i++) {
            const k = a.keys[i];
            if (k !== u) generator.add(k, aSets.get(k))
        }
        return generator.getSet();
    } else {
        const generator = Generator();
        for (let i = 0, _i = a.keys.length; i < _i; i++) {
            const k = a.keys[i];
            if (k === u) generator.add(k, OrderedSet.subtract(set, OrderedSet.ofSingleton(v)))
            else generator.add(k, aSets.get(k))
        }
        return generator.getSet();
    }
}

function subtractEE(a: AtomSetElements, b: AtomSetElements) {
    if (a === b) return Empty;

    const keysA = a.keys, keysB = b.keys;
    if (!SortedArray.areIntersecting(keysA, keysB)) return a;
    const r = SortedArray.findRange(keysA, SortedArray.min(keysB), SortedArray.max(keysB));
    const start = Interval.start(r), end = Interval.end(r);

    const generator = Generator();
    const { sets: aSets } = a;
    const { sets: bSets } = b;
    for (let i = 0; i < start; i++) {
        const k = keysA[i];
        generator.add(k, aSets.get(k));
    }
    for (let i = start; i < end; i++) {
        const k = keysA[i];
        const ak = aSets.get(k), bk = bSets.get(k);
        if (!!bk) {
            const subtraction = OrderedSet.subtract(ak, bk);
            generator.add(k, subtraction);
        } else {
            generator.add(k, ak);
        }
    }
    for (let i = end, _i = keysA.length; i < _i; i++) {
        const k = keysA[i];
        generator.add(k, aSets.get(k));
    }
    return generator.getSet();
}

function findUnion(sets: ArrayLike<AtomSetImpl>) {
    if (!sets.length) return Empty;
    if (sets.length === 1) return sets[0];
    if (sets.length === 2 && areEqual(sets[0], sets[1])) return sets[0];

    const eCount = { count: 0 };
    const ns = unionN(sets, eCount);
    if (!eCount.count) return ns;
    const generator = Generator();
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s !== 'number') unionInto(generator, s as AtomSetElements);
    }
    if (size(ns as AtomSetImpl) > 0) {
        if (typeof ns === 'number') unionIntoN(generator, ns as any);
        else unionInto(generator, ns as AtomSetElements);
    }
    return generator.getSet();
}

function unionN(sets: ArrayLike<AtomSetImpl>, eCount: { count: number }) {
    let countN = 0, countE = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        if (typeof sets[i] === 'number') countN++;
        else countE++;
    }
    eCount.count = countE;
    if (!countN) return Empty;
    if (countN === sets.length) return ofAtoms(sets as ArrayLike<Atom>);
    const packed = Atom.createEmptyArray(countN);
    let offset = 0;
    for (let i = 0, _i = sets.length; i < _i; i++) {
        const s = sets[i];
        if (typeof s === 'number') packed[offset++] = s;
    }
    return ofAtoms(packed as any);
}

function unionInto(builder: AtomSetGenerator, a: AtomSetElements) {
    const keys = a.keys;
    const { sets: aSets } = a;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        builder.addUnion(k, aSets.get(k));
    }
}

function unionIntoN(builder: AtomSetGenerator, a: Atom) {
    const u = Atom.unit(a);
    builder.addUnion(u, OrderedSet.ofSingleton(Atom.index(a)));
}