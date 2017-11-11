/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray, Interval, Iterator, OrderedSet as OS, IntMap } from 'mol-data/int'
import { sortArray } from 'mol-data/util/sort'
import { hash1 } from 'mol-data/util/hash-functions'
import Atom from '../../atom'
import AtomGroup from '../group'

/** Long and painful implementation starts here */

export type AtomSetImpl = { groups: IntMap<AtomGroup>, offsets: Int32Array, hashCode: number, keys: SortedArray }

export const Empty: AtomSetImpl = { groups: IntMap.Empty, offsets: new Int32Array(1), hashCode: 0, keys: SortedArray.Empty };

export function ofAtoms(atoms: ArrayLike<Atom>, template: AtomSetImpl): AtomSetImpl {
    return ofAtomsImpl(atoms, template);
}

export function getKeys(set: AtomSetImpl): SortedArray {
    return set.keys;
}

export function keyCount(set: AtomSetImpl): number {
    return set.keys.length;
}

export function hasKey(set: AtomSetImpl, key: number): boolean {
    return set.groups.has(key);
}

export function getKey(set: AtomSetImpl, index: number): number {
    return set.keys[index];
}

export function hasAtom(set: AtomSetImpl, t: Atom): boolean {
    const os = set.groups.get(Atom.unit(t));
    return !!os && AtomGroup.has(os, Atom.index(t));
}

export function getByKey(set: AtomSetImpl, key: number): AtomGroup {
    return set.groups.get(key) || AtomGroup.Empty;
}

export function getByIndex(set: AtomSetImpl, index: number): AtomGroup {
    const key = set.keys[index];
    return set.groups.get(key) || AtomGroup.Empty;
}

export function getAt(set: AtomSetImpl, i: number): Atom {
    const { offsets, keys } = set;
    const o = getOffsetIndex(offsets, i);
    if (o >= offsets.length - 1) return Atom.Zero;
    const k = keys[o];
    const e = AtomGroup.getAt(set.groups.get(k), i - offsets[o]);
    return Atom.create(k, e);
}

export function indexOf(set: AtomSetImpl, t: Atom) {
    const { keys } = set;
    const u = Atom.unit(t);
    const setIdx = SortedArray.indexOf(keys, u);
    if (setIdx < 0) return -1;
    const o = AtomGroup.indexOf(set.groups.get(u), Atom.index(t));
    if (o < 0) return -1;
    return set.offsets[setIdx] + o;
}

/** Number elements in the "child" sets */
export function size(set: AtomSetImpl) {
    return set.offsets[set.offsets.length - 1];
}

export function hashCode(set: AtomSetImpl) {
    if (set.hashCode !== -1) return set.hashCode;
    return computeHash(set);
}

export function areEqual(a: AtomSetImpl, b: AtomSetImpl) {
    if (a === b) return true;
    if (size(a) !== size(a)) return false;

    const keys = a.keys;
    if (!SortedArray.areEqual(keys, b.keys)) return false;
    const { groups: aG } = a;
    const { groups: bG } = b;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        if (!AtomGroup.areEqual(aG.get(k), bG.get(k))) return false;
    }
    return true;
}

export function areIntersecting(a: AtomSetImpl, b: AtomSetImpl) {
    if (a === b) return true;
    const keysA = a.keys, keysB = b.keys;
    if (!SortedArray.areIntersecting(a.keys, b.keys)) return false;
    const r = SortedArray.findRange(keysA, SortedArray.min(keysB), SortedArray.max(keysB));
    const start = Interval.start(r), end = Interval.end(r);
    const { groups: aG } = a;
    const { groups: bG } = b;
    for (let i = start; i < end; i++) {
        const k = keysA[i];
        const ak = aG.get(k), bk = bG.get(k);
        if (!!ak && !!bk && OS.areIntersecting(ak.atoms, bk.atoms)) return true;
    }
    return false;
}

export function intersect(a: AtomSetImpl, b: AtomSetImpl) {
    if (a === b) return a;

    const keysA = a.keys, keysB = b.keys;
    if (!SortedArray.areIntersecting(a.keys, b.keys)) return Empty;
    const r = SortedArray.findRange(keysA, SortedArray.min(keysB), SortedArray.max(keysB));
    const start = Interval.start(r), end = Interval.end(r);

    const { groups: aG } = a;
    const { groups: bG } = b;
    const generator = new ChildGenerator(a, b);
    for (let i = start; i < end; i++) {
        const k = keysA[i];
        const bk = bG.get(k);
        if (!bk) continue;
        const ak = aG.get(k);
        generator.add(k, AtomGroup.intersect(aG.get(k), bk), ak, bk);
    }
    return generator.getSet();
}

export function subtract(a: AtomSetImpl, b: AtomSetImpl) {
    if (a === b) return Empty;

    const keysA = a.keys, keysB = b.keys;
    if (!SortedArray.areIntersecting(keysA, keysB)) return a;
    const r = SortedArray.findRange(keysA, SortedArray.min(keysB), SortedArray.max(keysB));
    const start = Interval.start(r), end = Interval.end(r);

    const generator = new ChildGenerator(a, b);
    const { groups: aG } = a;
    const { groups: bG } = b;
    for (let i = 0; i < start; i++) {
        const k = keysA[i];
        const ak = aG.get(k);
        generator.addA(k, ak, ak);
    }
    for (let i = start; i < end; i++) {
        const k = keysA[i];
        const ak = aG.get(k), bk = bG.get(k);
        if (!!bk) {
            const subtraction = AtomGroup.subtract(ak, bk);
            generator.addA(k, subtraction, ak);
        } else {
            generator.addA(k, ak, ak);
        }
    }
    for (let i = end, _i = keysA.length; i < _i; i++) {
        const k = keysA[i];
        const ak = aG.get(k);
        generator.addA(k, ak, ak);
    }
    return generator.getSet();
}

export function unionMany(sets: ArrayLike<AtomSetImpl>, template: AtomSetImpl) {
    return findUnion(sets, template);
}

class ElementsIterator implements Iterator<Atom> {
    private unit: number = 0;
    private keyCount: number;
    private setIndex = -1;
    private currentIndex = 0;
    private currentSize = 0;
    private currentSet: OS = OS.Empty;

    hasNext: boolean = false;

    move() {
        if (!this.hasNext) return Atom.Zero;
        const ret = Atom.create(this.unit, OS.getAt(this.currentSet, this.currentIndex++));
        if (this.currentIndex >= this.currentSize) this.advance();
        return ret;
    }

    private advance() {
        if (++this.setIndex >= this.keyCount) {
            this.hasNext = false;
            return false;
        }
        this.unit = this.elements.keys[this.setIndex];
        this.currentSet = this.elements.groups.get(this.unit).atoms;
        this.currentIndex = 0;
        this.currentSize = OS.size(this.currentSet);
        return true;
    }

    constructor(private elements: AtomSetImpl) {
        this.keyCount = elements.keys.length;
        this.hasNext = this.keyCount > 0;
        this.advance();
    }
}

export function values(set: AtomSetImpl): Iterator<Atom> {
    return new ElementsIterator(set);
}

export class TemplateAtomSetGenerator {
    private keys: number[] = [];
    private groups = IntMap.Mutable<AtomGroup>();
    private templateGroups: IntMap<AtomGroup>;
    private equalGroups = 0;

    add(unit: number, group: AtomGroup) {
        if (AtomGroup.size(group) === 0) return;
        this.keys[this.keys.length] = unit;
        const templ = this.templateGroups.get(unit);
        if (AtomGroup.areEqual(templ, group)) {
            this.groups.set(unit, templ);
            this.equalGroups++;
        } else {
            this.groups.set(unit, group);
        }

    }

    getSet(): AtomSetImpl {
        if (this.equalGroups === this.template.keys.length) {
            return this.template;
        }
        return create(this.keys, this.groups);
    }

    constructor(private template: AtomSetImpl) {
        this.templateGroups = template.groups;
    }
}

export function TemplateGenerator(template: AtomSetImpl) {
    return new TemplateAtomSetGenerator(template);
}

export class AtomSetGenerator {
    private keys: number[] = [];
    private groups = IntMap.Mutable<AtomGroup>();

    add(unit: number, group: AtomGroup) {
        if (AtomGroup.size(group) === 0) return;
        this.keys[this.keys.length] = unit;
        this.groups.set(unit, group);
    }

    getSet(): AtomSetImpl {
        return create(this.keys, this.groups);
    }
}

export function Generator() {
    return new AtomSetGenerator();
}

/** When adding groups, compare them to existing ones. If they all match, return the whole original set. */
class ChildGenerator {
    private keys: number[] = [];
    private groups = IntMap.Mutable<AtomGroup>();
    private aEqual = 0;
    private bEqual = 0;

    add(unit: number, group: AtomGroup, a: AtomGroup, b: AtomGroup) {
        if (AtomGroup.size(group) === 0) return;
        if (a === group) this.aEqual++;
        if (b === group) this.bEqual++;
        this.keys[this.keys.length] = unit;
        this.groups.set(unit, group);
    }

    addA(unit: number, group: AtomGroup, a: AtomGroup) {
        if (AtomGroup.size(group) === 0) return;

        if (a === group) this.aEqual++;
        this.keys[this.keys.length] = unit;
        this.groups.set(unit, group);
    }

    constructor(private a: AtomSetImpl, private b: AtomSetImpl) {
    }

    getSet(): AtomSetImpl {
        if (this.aEqual === this.a.keys.length) return this.a;
        if (this.bEqual === this.b.keys.length) return this.b;
        return create(this.keys, this.groups);
    }
}

function create(keys: number[], groups: IntMap<AtomGroup>): AtomSetImpl {
    sortArray(keys);
    let runningSize = 0;
    const offsets = new Int32Array(keys.length + 1);
    for (let i = 0, _i = keys.length; i < _i; i++) {
        runningSize += AtomGroup.size(groups.get(keys[i]));
        offsets[i + 1] = runningSize;
    }
    return { keys: SortedArray.ofSortedArray(keys), groups: IntMap.asImmutable(groups), offsets, hashCode: -1 };
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

function ofAtomsImpl(xs: ArrayLike<Atom>, template: AtomSetImpl) {
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

    const generator = TemplateGenerator(template);
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        const group = AtomGroup.createNew(OS.ofSortedArray(normalizeArray(elements.get(k))));
        generator.add(k, group);
    }

    return generator.getSet();
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

function computeHash(set: AtomSetImpl) {
    const { keys, groups } = set;
    let hash = 23;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        hash = (31 * hash + k) | 0;
        hash = (31 * hash + AtomGroup.hashCode(groups.get(k))) | 0;
    }
    hash = (31 * hash + size(set)) | 0;
    hash = hash1(hash);
    set.hashCode = hash;
    return hash;
}

function findUnion(sets: ArrayLike<AtomSetImpl>, template: AtomSetImpl) {
    if (!sets.length) return Empty;
    if (sets.length === 1) return sets[0];
    if (sets.length === 2 && sets[0] === sets[1]) return sets[0];

    const keys: number[] = [];
    const groups = IntMap.Mutable<AtomGroup>();
    for (let i = 0, _i = sets.length; i < _i; i++) {
        unionInto(keys, groups, sets[i]);
    }

    return normalizeUnion(keys, groups, template);
}

function normalizeUnion(keys: number[], groups: IntMap.Mutable<AtomGroup>, template: AtomSetImpl) {
    let equalCount = 0;
    let tg = template.groups;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        const a = groups.get(k), t = tg.get(k);
        if (AtomGroup.areEqual(a, t)) {
            groups.set(k, t);
            equalCount++;
        }
    }
    return equalCount === template.keys.length ? template : create(keys, groups);
}

function unionInto(keys: number[], groups: IntMap.Mutable<AtomGroup>, a: AtomSetImpl) {
    const setKeys = a.keys;
    const { groups: aG } = a;
    for (let i = 0, _i = setKeys.length; i < _i; i++) {
        const k = setKeys[i];
        if (groups.has(k)) {
            groups.set(k, AtomGroup.union(aG.get(k), groups.get(k)))
        } else {
            groups.set(k, aG.get(k));
        }
    }
}