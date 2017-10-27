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
 * while (it.hasNext) { const v = it.move(); ... } 
 */
interface Iterator<T> {
    readonly hasNext: boolean,
    move(): T
}

class ArrayIteratorImpl<T> implements Iterator<T> {
    private xs: ArrayLike<T> = [];
    private index: number = -1;
    private length: number = 0;
    private lastValue: T;

    hasNext: boolean;

    move() {
        ++this.index;
        this.lastValue = this.xs[this.index];
        this.hasNext = this.index < this.length - 1;
        return this.lastValue;
    }

    constructor(xs: ArrayLike<T>) {
        this.length = xs.length;
        this.hasNext = xs.length > 0;
        this.xs = xs;
        this.index = -1;
        // try to avoid deoptimization with undefined values
        this.lastValue = xs.length > 0 ? xs[0] : void 0 as any;
    }
}

class RangeIteratorImpl implements Iterator<number> {
    private value: number;

    hasNext: boolean;

    move() {
        ++this.value;
        this.hasNext = this.value < this.max;
        return this.value;
    }

    constructor(min: number, private max: number) {
        this.value = min - 1;
        this.hasNext = max >= min;
    }
}

class ValueIterator<T> implements Iterator<T> {
    hasNext = true;
    move() { this.hasNext = false; return this.value; }
    constructor(private value: T) { }
}

namespace Iterator {
    export const Empty: Iterator<any> = new RangeIteratorImpl(0, -1);
    export function Array<T>(xs: ArrayLike<T>): Iterator<T> { return new ArrayIteratorImpl<T>(xs); }
    export function Value<T>(value: T): Iterator<T> { return new ValueIterator(value); }
    export function Range(min: number, max: number): Iterator<number> { return new RangeIteratorImpl(min, max); }
}

export default Iterator