/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray, Interval, Iterator, OrderedSet as OS, IntMap } from 'mol-data/int'
import { sortArray } from 'mol-data/util/sort'
import { hash1 } from 'mol-data/util/hash-functions'
import Element from '../../element'
import ElementGroup from '../group'

/** Long and painful implementation starts here */

export type ElementSetImpl = { groups: IntMap<ElementGroup>, offsets: Int32Array, hashCode: number, keys: SortedArray }

export const Empty: ElementSetImpl = { groups: IntMap.Empty, offsets: new Int32Array(1), hashCode: 0, keys: SortedArray.Empty };

export function ofElements(elements: ArrayLike<Element>, template: ElementSetImpl): ElementSetImpl {
    return ofElementsImpl(elements, template);
}

export function singleton(element: Element, template: ElementSetImpl) {
    return singletonImpl(element, template);
}

export function getKeys(set: ElementSetImpl): SortedArray {
    return set.keys;
}

export function keyCount(set: ElementSetImpl): number {
    return set.keys.length;
}

export function hasKey(set: ElementSetImpl, key: number): boolean {
    return set.groups.has(key);
}

export function getKey(set: ElementSetImpl, index: number): number {
    return set.keys[index];
}

export function hasAtom(set: ElementSetImpl, t: Element): boolean {
    const os = set.groups.get(Element.unit(t));
    return !!os && ElementGroup.has(os, Element.index(t));
}

export function getByKey(set: ElementSetImpl, key: number): ElementGroup {
    return set.groups.get(key) || ElementGroup.Empty;
}

export function getByIndex(set: ElementSetImpl, index: number): ElementGroup {
    const key = set.keys[index];
    return set.groups.get(key) || ElementGroup.Empty;
}

export function getAt(set: ElementSetImpl, i: number): Element {
    const { offsets, keys } = set;
    const o = getOffsetIndex(offsets, i);
    if (o >= offsets.length - 1) return Element.Zero;
    const k = keys[o];
    const e = ElementGroup.getAt(set.groups.get(k), i - offsets[o]);
    return Element.create(k, e);
}

export function indexOf(set: ElementSetImpl, t: Element) {
    const { keys } = set;
    const u = Element.unit(t);
    const setIdx = SortedArray.indexOf(keys, u);
    if (setIdx < 0) return -1;
    const o = ElementGroup.indexOf(set.groups.get(u), Element.index(t));
    if (o < 0) return -1;
    return set.offsets[setIdx] + o;
}

/** Number elements in the "child" sets */
export function size(set: ElementSetImpl) {
    return set.offsets[set.offsets.length - 1];
}

export function hashCode(set: ElementSetImpl) {
    if (set.hashCode !== -1) return set.hashCode;
    return computeHash(set);
}

export function areEqual(a: ElementSetImpl, b: ElementSetImpl) {
    if (a === b) return true;
    if (size(a) !== size(a)) return false;

    const keys = a.keys;
    if (!SortedArray.areEqual(keys, b.keys)) return false;
    const { groups: aG } = a;
    const { groups: bG } = b;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        if (!ElementGroup.areEqual(aG.get(k), bG.get(k))) return false;
    }
    return true;
}

export function areIntersecting(a: ElementSetImpl, b: ElementSetImpl) {
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
        if (!!ak && !!bk && OS.areIntersecting(ak.elements, bk.elements)) return true;
    }
    return false;
}

export function intersect(a: ElementSetImpl, b: ElementSetImpl) {
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
        generator.add(k, ElementGroup.intersect(aG.get(k), bk), ak, bk);
    }
    return generator.getSet();
}

export function subtract(a: ElementSetImpl, b: ElementSetImpl) {
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
            const subtraction = ElementGroup.subtract(ak, bk);
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

export function unionMany(sets: ArrayLike<ElementSetImpl>, template: ElementSetImpl) {
    return findUnion(sets, template);
}

class ElementsIterator implements Iterator<Element> {
    private unit: number = 0;
    private keyCount: number;
    private setIndex = -1;
    private currentIndex = 0;
    private currentSize = 0;
    private currentSet: OS = OS.Empty;

    hasNext: boolean = false;

    move() {
        if (!this.hasNext) return Element.Zero;
        const ret = Element.create(this.unit, OS.getAt(this.currentSet, this.currentIndex++));
        if (this.currentIndex >= this.currentSize) this.advance();
        return ret;
    }

    private advance() {
        if (++this.setIndex >= this.keyCount) {
            this.hasNext = false;
            return false;
        }
        this.unit = this.elements.keys[this.setIndex];
        this.currentSet = this.elements.groups.get(this.unit).elements;
        this.currentIndex = 0;
        this.currentSize = OS.size(this.currentSet);
        return true;
    }

    constructor(private elements: ElementSetImpl) {
        this.keyCount = elements.keys.length;
        this.hasNext = this.keyCount > 0;
        this.advance();
    }
}

export function values(set: ElementSetImpl): Iterator<Element> {
    return new ElementsIterator(set);
}

export class TemplateAtomSetGenerator {
    private keys: number[] = [];
    private groups = IntMap.Mutable<ElementGroup>();
    private templateGroups: IntMap<ElementGroup>;
    private equalGroups = 0;

    add(unit: number, set: OS) {
        if (OS.size(set) === 0) return;
        this.keys[this.keys.length] = unit;
        if (this.templateGroups.has(unit)) {
            const t = this.templateGroups.get(unit);
            if (OS.areEqual(t.elements, set)) {
                this.groups.set(unit, t);
                this.equalGroups++;
            } else {
                this.groups.set(unit, ElementGroup.createNew(set));
            }
        } else {
            this.groups.set(unit, ElementGroup.createNew(set));
        }
    }

    getSet(): ElementSetImpl {
        if (this.equalGroups === this.template.keys.length && this.equalGroups === this.keys.length) {
            return this.template;
        }
        return create(this.keys, this.groups);
    }

    constructor(private template: ElementSetImpl) {
        this.templateGroups = template.groups;
    }
}

export function TemplateGenerator(template: ElementSetImpl) {
    return new TemplateAtomSetGenerator(template);
}

export class AtomSetGenerator {
    private keys: number[] = [];
    private groups = IntMap.Mutable<ElementGroup>();

    add(unit: number, group: ElementGroup) {
        if (ElementGroup.size(group) === 0) return;
        this.keys[this.keys.length] = unit;
        this.groups.set(unit, group);
    }

    getSet(): ElementSetImpl {
        return create(this.keys, this.groups);
    }
}

export function Generator() {
    return new AtomSetGenerator();
}

/** When adding groups, compare them to existing ones. If they all match, return the whole original set. */
class ChildGenerator {
    private keys: number[] = [];
    private groups = IntMap.Mutable<ElementGroup>();
    private aEqual = 0;
    private bEqual = 0;

    add(unit: number, group: ElementGroup, a: ElementGroup, b: ElementGroup) {
        if (ElementGroup.size(group) === 0) return;
        if (a === group) this.aEqual++;
        if (b === group) this.bEqual++;
        this.keys[this.keys.length] = unit;
        this.groups.set(unit, group);
    }

    addA(unit: number, group: ElementGroup, a: ElementGroup) {
        if (ElementGroup.size(group) === 0) return;

        if (a === group) this.aEqual++;
        this.keys[this.keys.length] = unit;
        this.groups.set(unit, group);
    }

    constructor(private a: ElementSetImpl, private b: ElementSetImpl) {
    }

    getSet(): ElementSetImpl {
        if (this.aEqual === this.a.keys.length) return this.a;
        if (this.bEqual === this.b.keys.length) return this.b;
        return create(this.keys, this.groups);
    }
}

function create(keys: number[], groups: IntMap<ElementGroup>): ElementSetImpl {
    sortArray(keys);
    let runningSize = 0;
    const offsets = new Int32Array(keys.length + 1);
    for (let i = 0, _i = keys.length; i < _i; i++) {
        runningSize += ElementGroup.size(groups.get(keys[i]));
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

function ofElementsImpl(xs: ArrayLike<Element>, template: ElementSetImpl) {
    if (xs.length === 0) return Empty;

    const elements = IntMap.Mutable<number[]>();
    const keys: number[] = [];
    for (let i = 0, _i = xs.length; i < _i; i++) {
        const x = xs[i];
        const u = Element.unit(x), v = Element.index(x);
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
        generator.add(k, OS.ofSortedArray(normalizeArray(elements.get(k))));
    }

    return generator.getSet();
}

function singletonImpl(element: Element, template: ElementSetImpl) {
    const k = Element.unit(element), i = Element.index(element);
    const { groups } = template;
    const gs = IntMap.Mutable<ElementGroup>();
    if (groups.has(k)) {
        const g = groups.get(k);
        if (ElementGroup.size(g) === 1 && ElementGroup.getAt(g, 0) === i) {
            gs.set(k, g);
            return create([k], gs);
        }
    }
    gs.set(k, ElementGroup.createNew(OS.ofSingleton(i)));
    return create([k], gs);
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

function computeHash(set: ElementSetImpl) {
    const { keys, groups } = set;
    let hash = 23;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        hash = (31 * hash + k) | 0;
        hash = (31 * hash + ElementGroup.hashCode(groups.get(k))) | 0;
    }
    hash = (31 * hash + size(set)) | 0;
    hash = hash1(hash);
    set.hashCode = hash;
    return hash;
}

function findUnion(sets: ArrayLike<ElementSetImpl>, template: ElementSetImpl) {
    if (!sets.length) return Empty;
    if (sets.length === 1) return sets[0];
    if (sets.length === 2 && sets[0] === sets[1]) return sets[0];

    const keys: number[] = [];
    const groups = IntMap.Mutable<ElementGroup>();
    for (let i = 0, _i = sets.length; i < _i; i++) {
        unionInto(keys, groups, sets[i]);
    }

    return normalizeUnion(keys, groups, template);
}

function normalizeUnion(keys: number[], groups: IntMap.Mutable<ElementGroup>, template: ElementSetImpl) {
    let equalCount = 0;
    let tg = template.groups, a: ElementGroup, t: ElementGroup;
    for (let i = 0, _i = keys.length; i < _i; i++) {
        const k = keys[i];
        if (tg.has(k) && ElementGroup.areEqual(a = groups.get(k), t = tg.get(k))) {
            groups.set(k, t);
            equalCount++;
        }
    }
    return equalCount === template.keys.length && equalCount === keys.length ? template : create(keys, groups);
}

function unionInto(keys: number[], groups: IntMap.Mutable<ElementGroup>, a: ElementSetImpl) {
    const setKeys = a.keys;
    const { groups: aG } = a;
    for (let i = 0, _i = setKeys.length; i < _i; i++) {
        const k = setKeys[i];
        if (groups.has(k)) {
            groups.set(k, ElementGroup.union(aG.get(k), groups.get(k)))
        } else {
            keys[keys.length] = k;
            groups.set(k, aG.get(k));
        }
    }
}