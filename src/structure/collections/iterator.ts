/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * ES6 compatible mutable iterator. Use with care. 
 *
 * "Idiomatic" usage:
 *
 * const it = ...;
 * for (let v = it.nextValue(); !it.done; v = it.nextValue()) { ... }
 */
interface Iterator<T> {
    [Symbol.iterator](): Iterator<T>,
    readonly done: boolean,
    readonly value: T,
    next(): { done: boolean, value: T },
    nextValue(): T
}

class ArrayIteratorImpl<T> implements Iterator<T> {
    private xs: ArrayLike<T> = [];
    private index: number = -1;
    private length: number = 0;

    [Symbol.iterator]() { return this; };
    done = true;
    value: T = void 0 as any;

    next() {
        const index = ++this.index;
        if (index < this.length) this.value = this.xs[index];
        else this.done = true;
        return this;
    }

    nextValue() { return this.next().value; }

    constructor(xs: ArrayLike<T>) {
        this.length = xs.length;
        this.done = false;
        this.xs = xs;
        this.index = -1;
        return this;
    }
}

class RangeIteratorImpl implements Iterator<number> {
    [Symbol.iterator]() { return this; };
    done = true;
    value: number;

    next() {
        ++this.value;
        this.done = this.value > this.max;
        return this;
    }

    nextValue() { return this.next().value;  }

    constructor(min: number, private max: number) {
        this.value = min - 1;
        this.done = false;
        return this;
    }
}

namespace Iterator {
    export const Empty: Iterator<any> = new RangeIteratorImpl(0, -1);
    export function Array<T>(xs: ArrayLike<T>): Iterator<T> { return new ArrayIteratorImpl<T>(xs); }
    export function Value(value: number): Iterator<number> { return new RangeIteratorImpl(value, value); }
    export function Range(min: number, max: number): Iterator<number> { return new RangeIteratorImpl(min, max); }

    export function toArray<T>(it: Iterator<T>): T[] {
        const ret = [];
        for (let v = it.nextValue(); !it.done; v = it.nextValue()) ret[ret.length] = v;
        return ret;
    }
}

export default Iterator