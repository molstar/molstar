/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// adpated from https://github.com/dankogai/js-combinatorics, MIT 2013-2016 Dan Kogai

import Iterator from '../iterator';

function P(m: number, n: number) {
    let p = 1;
    while (n--) p *= m--;
    return p;
}

function C(m: number, n: number) {
    if (n > m) return 0;
    return P(m, n) / P(n, n);
}

function nextIndex(n: number) {
    const smallest = n & -n;
    const ripple = n + smallest;
    const newSmallest = ripple & -ripple;
    const ones = ((newSmallest / smallest) >> 1) - 1;
    return ripple | ones;
};

export class CombinationIterator<T> implements Iterator<ReadonlyArray<T>> {
    private value: T[]
    private index: number
    private maxIndex: number

    size: number
    hasNext: boolean = false;

    move() {
        if (this.hasNext) {
            let i = 0, j = 0, n = this.index;
            for (; n; n >>>= 1, i++) {
                if (n & 1) this.value[j++] = this.array[i];
            }
            this.index = nextIndex(this.index);
            this.hasNext = this.index < this.maxIndex;
        }
        return this.value;
    }

    constructor(private array: T[], count: number) {
        this.index = (1 << count) - 1;
        this.size = C(array.length, count);
        this.maxIndex = 1 << array.length,

        this.value = new Array(count);
        this.hasNext = count > 0 && count <= array.length;
    }
}

export function combinations<T>(array: T[], count: number): T[][] {
    const out: T[][] = [];
    const combinationIt = new CombinationIterator(array, count);
    while (combinationIt.hasNext) out.push(combinationIt.move().slice());
    return out;
}