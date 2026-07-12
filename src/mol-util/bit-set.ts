/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

/**
 * A dense, fixed-size bit set backed by a Uint32Array. Bits are addressed 0..nbits-1.
 *
 * Intentionally minimal — extend as needed. The in-place set operations (and/andNot/or) assume the
 * operand has the same capacity.
 */
export class BitSet {
    readonly nbits: number;
    readonly words: Uint32Array;

    constructor(nbits: number) {
        this.nbits = nbits;
        this.words = new Uint32Array((nbits + 31) >>> 5);
    }

    /** A set with every bit in [0, nbits) set. */
    static full(nbits: number): BitSet {
        const b = new BitSet(nbits);
        b.words.fill(0xffffffff);
        const rem = nbits & 31;
        if (rem !== 0) b.words[b.words.length - 1] = ~(0xffffffff << rem);
        return b;
    }

    clone(): BitSet {
        const b = new BitSet(this.nbits);
        b.words.set(this.words);
        return b;
    }

    /** Clear all bits. */
    zero() { this.words.fill(0); }

    set(i: number) { this.words[i >>> 5] |= (1 << (i & 31)); }
    clear(i: number) { this.words[i >>> 5] &= ~(1 << (i & 31)); }
    get(i: number): boolean { return (this.words[i >>> 5] & (1 << (i & 31))) !== 0; }

    /** In-place intersection with `o` (same capacity). */
    and(o: BitSet) { const w = this.words, x = o.words; for (let k = 0, n = w.length; k < n; k++) w[k] &= x[k]; }
    /** In-place set difference with `o` (same capacity). */
    andNot(o: BitSet) { const w = this.words, x = o.words; for (let k = 0, n = w.length; k < n; k++) w[k] &= ~x[k]; }
    /** In-place union with `o` (same capacity). */
    or(o: BitSet) { const w = this.words, x = o.words; for (let k = 0, n = w.length; k < n; k++) w[k] |= x[k]; }

    isEmpty(): boolean {
        const w = this.words;
        for (let k = 0, n = w.length; k < n; k++) if (w[k] !== 0) return false;
        return true;
    }

    /** Number of set bits. */
    cardinality(): number {
        const w = this.words;
        let c = 0;
        for (let k = 0, n = w.length; k < n; k++) c += bitCount(w[k]);
        return c;
    }

    /** Index of the first set bit at or after `from`, or -1 if there is none. */
    nextSetBit(from: number): number {
        if (from < 0) from = 0;
        if (from >= this.nbits) return -1;
        const w = this.words;
        let wi = from >>> 5;
        let word = w[wi] & (0xffffffff << (from & 31));
        while (true) {
            if (word !== 0) {
                const bit = (wi << 5) + (31 - Math.clz32(word & -word));
                return bit < this.nbits ? bit : -1;
            }
            if (++wi >= w.length) return -1;
            word = w[wi];
        }
    }
}

/** Population count (number of set bits) of a 32-bit word. */
function bitCount(n: number): number {
    n = n - ((n >>> 1) & 0x55555555);
    n = (n & 0x33333333) + ((n >>> 2) & 0x33333333);
    return (((n + (n >>> 4)) & 0x0f0f0f0f) * 0x01010101) >>> 24;
}