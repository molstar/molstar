/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * ES6 compatible iterator.
 *
 * "Idiomatic" usage is to use the move function, because it does not make any allocations.
 *
 * const it = ...;
 * for (let v = it.move(); it.hasNext; v = it.move()) { ... }
 */
interface Iterator<T> {
    [Symbol.iterator](): Iterator<T>,
    readonly hasNext: boolean,
    next(): { done: boolean, value: T },
    move(): T
}

class ArrayIteratorImpl<T> implements Iterator<T> {
    private xs: ArrayLike<T> = [];
    private index: number = -1;
    private length: number = 0;
    private lastValue: T;

    [Symbol.iterator]() { return this; };
    hasNext: boolean;

    next() {
        const value = this.move();
        return { value, done: !this.hasNext };
    }

    move() {
        const index = ++this.index;
        if (index < this.length) this.lastValue = this.xs[index];
        else this.hasNext = false;
        return this.lastValue;
    }

    constructor(xs: ArrayLike<T>) {
        this.length = xs.length;
        this.hasNext = xs.length > 0;
        this.xs = xs;
        this.index = -1;
        // try to avoid deoptimization with undefined values
        this.lastValue = xs.length > 0 ? xs[0] : void 0 as any;
        return this;
    }
}

class RangeIteratorImpl implements Iterator<number> {
    private value: number;

    [Symbol.iterator]() { return this; };
    hasNext: boolean;

    next() {
        const value = this.value;
        return { value, done: !this.hasNext }
    }

    move() {
        ++this.value;
        this.hasNext = this.value <= this.max;
        return this.value;
    }

    constructor(min: number, private max: number) {
        this.value = min - 1;
        this.hasNext = max >= min;
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
        for (let v = it.move(); it.hasNext; v = it.move()) ret[ret.length] = v;
        return ret;
    }
}

export default Iterator