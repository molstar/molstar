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
 * for (let v = it.reset(data).nextValue(); !it.done; v = it.nextValue()) { ... }
 */
interface Iterator<T, Data = any> {
    [Symbol.iterator](): Iterator<T, Data>,
    readonly done: boolean,
    readonly value: T,
    next(): { done: boolean, value: T },

    reset(data: Data): Iterator<T, Data>,
    nextValue(): T
}

class __EmptyIterator implements Iterator<any, any> { // tslint:disable-line:class-name
    [Symbol.iterator]() { return this; }
    done = true;
    value = void 0;
    next() { return this; }
    nextValue() { return this.value; }
    reset(value: undefined) { return this; }
}

class __SingletonIterator<T> implements Iterator<T, T> { // tslint:disable-line:class-name
    private yielded = false;

    [Symbol.iterator]() { return this; }
    done = false;
    value: T;
    next() { this.done = this.yielded; this.yielded = true; return this; }
    nextValue() { return this.next().value; }
    reset(value: T) { this.value = value; this.done = false; this.yielded = false; return this; }

    constructor(value: T) { this.value = value; }
}


class __ArrayIterator<T> implements Iterator<T, ArrayLike<T>> { // tslint:disable-line:class-name
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

    reset(xs: ArrayLike<T>) {
        this.length = xs.length;
        this.done = false;
        this.xs = xs;
        this.index = -1;
        return this;
    }
}

type Range = { min: number, max: number }
class __RangeIterator implements Iterator<number, Range> { // tslint:disable-line:class-name
    private min: number;
    private max: number;

    [Symbol.iterator]() { return this; };
    done = true;
    value: number;

    next() {
        ++this.value;
        this.done = this.value >= this.max;
        return this;
    }

    nextValue() { return this.next().value;  }

    reset({ min, max}: Range) {
        this.min = min;
        this.max = max;
        this.value = min - 1;
        this.done = false;
        return this;
    }

    constructor(bounds: Range) { this.reset(bounds); }
}

export const EmptyIterator: Iterator<any> = new __EmptyIterator();
export function SingletonIterator<T>(value: T): Iterator<T, T> { return new __SingletonIterator(value); }
export function ArrayIterator<T>(xs?: ArrayLike<T>): Iterator<T, ArrayLike<T>> {
    const ret = new __ArrayIterator<T>();
    if (xs) ret.reset(xs);
    return ret;
}
export function RangeIterator(bounds?: Range): Iterator<number, Range> { return new __RangeIterator(bounds || { min: 0, max: 0 }); }

export function toArray<T>(it: Iterator<T>): T[] {
    const ret = [];
    for (let v = it.nextValue(); !it.done; v = it.nextValue()) ret[ret.length] = v;
    return ret;
}

export default Iterator