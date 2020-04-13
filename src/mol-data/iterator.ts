/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * "Idiomatic" usage:
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

    hasNext: boolean = false;

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
    private value: number = 0;
    hasNext: boolean = false;

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

class MapIteratorImpl<T, R> implements Iterator<R> {
    hasNext: boolean = false;

    move() {
        const v = this.f(this.base.move());
        this.hasNext = this.base.hasNext;
        return v;
    }

    constructor(private base: Iterator<T>, private f: (v: T) => R) {
        this.hasNext = base.hasNext;
    }
}

class FilterIteratorImpl<T> implements Iterator<T> {
    private next: T;
    hasNext: boolean;

    move() {
        const ret = this.next;
        this.hasNext = this.findNext();
        return ret;
    }

    private findNext() {
        while (this.base.hasNext) {
            this.next = this.base.move();
            if (this.p(this.next)) return true;
        }
        return false;
    }

    constructor(private base: Iterator<T>, private p: (v: T) => boolean) {
        this.hasNext = this.findNext();
    }
}

namespace Iterator {
    export const Empty: Iterator<any> = new RangeIteratorImpl(0, -1);
    export function Array<T>(xs: ArrayLike<T>): Iterator<T> { return new ArrayIteratorImpl<T>(xs); }
    export function Value<T>(value: T): Iterator<T> { return new ValueIterator(value); }
    export function Range(min: number, max: number): Iterator<number> { return new RangeIteratorImpl(min, max); }
    export function map<T, R>(base: Iterator<T>, f: (v: T) => R): Iterator<R> { return new MapIteratorImpl(base, f); }
    export function filter<T>(base: Iterator<T>, p: (v: T) => boolean): Iterator<T> { return new FilterIteratorImpl(base, p); }

    // Iterate until first truthy value is returned.
    export function forEach<T, Ctx>(it: Iterator<T>, f: (v: T, ctx: Ctx) => any, ctx: Ctx): Ctx {
        while (it.hasNext) {
            const c = f(it.move(), ctx);
            if (c) return ctx;
        }
        return ctx;
    }
}

export default Iterator;