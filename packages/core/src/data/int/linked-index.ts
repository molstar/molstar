/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A data structure useful for graph traversal */
interface LinkedIndex {
    readonly head: number,
    has(i: number): boolean,
    remove(i: number): void
}

function LinkedIndex(size: number): LinkedIndex {
    return new LinkedIndexImpl(size);
}

class LinkedIndexImpl implements LinkedIndex {
    private prev: Int32Array;
    private next: Int32Array;
    head: number;

    remove(i: number) {
        const { prev, next } = this;
        const p = prev[i], n = next[i];
        if (p >= 0) {
            next[p] = n;
            prev[i] = -1;
        }
        if (n >= 0) {
            prev[n] = p;
            next[i] = -1;
        }
        if (i === this.head) {
            if (p < 0) this.head = n;
            else this.head = p;
        }
    }

    has(i: number) {
        return this.prev[i] >= 0 || this.next[i] >= 0 || this.head === i;
    }

    constructor(size: number) {
        this.head = size > 0 ? 0 : -1;
        this.prev = new Int32Array(size);
        this.next = new Int32Array(size);

        for (let i = 0; i < size; i++) {
            this.next[i] = i + 1;
            this.prev[i] = i - 1;
        }
        this.prev[0] = -1;
        this.next[size - 1] = -1;
    }
}

export { LinkedIndex };