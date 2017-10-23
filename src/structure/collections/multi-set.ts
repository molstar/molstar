/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from './ordered-set'
import Iterator from './iterator'
import IntPair from './int-pair'
import { sortArray } from './sort'

type MultiSetElements = { [id: number]: OrderedSet, keys: OrderedSet }
type MultiSet = number | MultiSetElements

namespace MultiSet {
    export const Empty: MultiSet = { keys: OrderedSet.Empty };

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

    function isArrayLike(x: any): x is ArrayLike<number> {
        return x && (typeof x.length === 'number' && (x instanceof Array || !!x.buffer));
    }

    export function create(data: number | ArrayLike<number> | IntPair | { [id: number]: OrderedSet }): MultiSet {
        if (typeof data === 'number') return data;
        if (IntPair.is(data)) return IntPair.pack(data);
        if (isArrayLike(data)) return ofPackedPairs(data);
        const keys = [];
        for (const _k of Object.keys(data)) {
            const k = +_k;
            if (data[k].size > 0) keys[keys.length] = k;
        }
        if (!keys.length) return Empty;
        sortArray(keys);
        const ret = Object.create(null);
        ret.keys = OrderedSet.ofSortedArray(keys);
        for (const k of keys) ret[k] = data[k];
        return ret;
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

    function ofPackedPairs(xs: ArrayLike<number>): MultiSet {
        if (xs.length === 0) return Empty;
        const sets: { [key: number]: number[] } = Object.create(null);
        const p = IntPair.zero();
        for (let i = 0, _i = xs.length; i < _i; i++) {
            IntPair.unpack(xs[i], p);
            const set = sets[p.fst];
            if (set) set[set.length] = p.snd;
            else sets[p.fst] = [p.snd];
        }
        const ret: { [key: number]: OrderedSet } = Object.create(null);
        for (const _k of Object.keys(sets)) {
            const k = +_k;
            ret[k] = OrderedSet.ofSortedArray(sortArray(sets[k]));
        }
        return create(ret);
    }

    // TODO: size, ofObject, hashCode (and cache it), equal, union, intersection, subtraction

    export function union(sets: ArrayLike<MultiSet>): MultiSet {
        return 0 as any;
    }
}

export default MultiSet