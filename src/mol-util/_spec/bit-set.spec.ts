/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Sebastian Bittrich <sebastian.m.bittrich@gmail.com>
 */

import { BitSet } from '../bit-set';

function setBits(b: BitSet): number[] {
    const out: number[] = [];
    for (let i = b.nextSetBit(0); i >= 0; i = b.nextSetBit(i + 1)) out.push(i);
    return out;
}

describe('BitSet', () => {
    it('set/get/clear including word boundaries and the sign bit', () => {
        const b = new BitSet(70);
        for (const i of [0, 31, 32, 63, 64, 69]) {
            expect(b.get(i)).toBe(false);
            b.set(i);
            expect(b.get(i)).toBe(true);
        }
        b.clear(32);
        expect(b.get(32)).toBe(false);
        expect(b.get(31)).toBe(true); // neighbour untouched
        expect(setBits(b)).toEqual([0, 31, 63, 64, 69]);
    });

    it('allocates ceil(nbits/32) words', () => {
        expect(new BitSet(0).words.length).toBe(0);
        expect(new BitSet(1).words.length).toBe(1);
        expect(new BitSet(32).words.length).toBe(1);
        expect(new BitSet(33).words.length).toBe(2);
    });

    it('cardinality / isEmpty', () => {
        const b = new BitSet(40);
        expect(b.isEmpty()).toBe(true);
        expect(b.cardinality()).toBe(0);
        b.set(5); b.set(39); b.set(5); // idempotent
        expect(b.isEmpty()).toBe(false);
        expect(b.cardinality()).toBe(2);
        b.zero();
        expect(b.isEmpty()).toBe(true);
        expect(b.cardinality()).toBe(0);
    });

    it('full() sets exactly [0, nbits) even when not a multiple of 32', () => {
        const b = new BitSet(40);
        expect(b.cardinality()).toBe(0);
        const f = BitSet.full(40);
        expect(f.cardinality()).toBe(40);
        expect(f.get(39)).toBe(true);
        expect(f.get(40)).toBe(false); // beyond nbits, within the last word
        expect(f.get(63)).toBe(false);
        expect(setBits(f)).toEqual(Array.from({ length: 40 }, (_, i) => i));

        const exact = BitSet.full(64);
        expect(exact.cardinality()).toBe(64);
    });

    it('clone is an independent copy', () => {
        const a = new BitSet(64);
        a.set(1); a.set(40);
        const c = a.clone();
        c.set(2);
        expect(setBits(a)).toEqual([1, 40]);
        expect(setBits(c)).toEqual([1, 2, 40]);
    });

    it('in-place and / andNot / or', () => {
        const mk = (...bits: number[]) => { const b = new BitSet(64); for (const i of bits) b.set(i); return b; };

        const a1 = mk(1, 2, 3, 40);
        a1.and(mk(2, 3, 4, 41));
        expect(setBits(a1)).toEqual([2, 3]);

        const a2 = mk(1, 2, 3, 40);
        a2.andNot(mk(2, 40));
        expect(setBits(a2)).toEqual([1, 3]);

        const a3 = mk(1, 40);
        a3.or(mk(2, 40, 63));
        expect(setBits(a3)).toEqual([1, 2, 40, 63]);
    });

    it('nextSetBit crosses words, handles the high bit, and terminates', () => {
        const b = new BitSet(100);
        expect(b.nextSetBit(0)).toBe(-1); // empty
        b.set(31); b.set(64); b.set(99);
        expect(b.nextSetBit(0)).toBe(31);
        expect(b.nextSetBit(31)).toBe(31);
        expect(b.nextSetBit(32)).toBe(64); // skips the rest of word 0 and all of word 1
        expect(b.nextSetBit(65)).toBe(99);
        expect(b.nextSetBit(100)).toBe(-1); // at/after nbits
        expect(b.nextSetBit(-5)).toBe(31); // negative clamps to 0
    });
});