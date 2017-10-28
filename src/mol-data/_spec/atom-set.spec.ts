/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from '../../mol-base/collections/integer/ordered-set'
import AtomSet from '../atom-set'
import Atom from '../atom'

describe('atom set', () => {
    const p = (i: number, j: number) => Atom.create(i, j);

    function setToPairs(set: AtomSet): ArrayLike<Atom> {
        const ret: Atom[] = [];
        const it = AtomSet.atoms(set);
        while (it.hasNext) {
            ret[ret.length] = it.move();
        }
        return ret;
    }

    it('singleton pair', () => {
        const set = AtomSet.create(p(10, 11));
        expect(setToPairs(set)).toEqual([p(10, 11)]);
        expect(AtomSet.atomHas(set, p(10, 11))).toBe(true);
        expect(AtomSet.atomHas(set, p(11, 11))).toBe(false);
        expect(AtomSet.atomGetAt(set, 0)).toBe(p(10, 11));
        expect(AtomSet.atomCount(set)).toBe(1);
    });

    it('singleton number', () => {
        const set = AtomSet.create(p(10, 11));
        expect(setToPairs(set)).toEqual([p(10, 11)]);
    });

    it('multi', () => {
        const set = AtomSet.create({
            1: OrderedSet.ofSortedArray([4, 6, 7]),
            3: OrderedSet.ofRange(0, 1),
        });
        const ret = [p(1, 4), p(1, 6), p(1, 7), p(3, 0), p(3, 1)];
        expect(AtomSet.atomCount(set)).toBe(ret.length);
        expect(setToPairs(set)).toEqual([p(1, 4), p(1, 6), p(1, 7), p(3, 0), p(3, 1)]);
        expect(AtomSet.atomHas(set, p(10, 11))).toBe(false);
        expect(AtomSet.atomHas(set, p(3, 0))).toBe(true);
        expect(AtomSet.atomHas(set, p(1, 7))).toBe(true);
        for (let i = 0; i < AtomSet.atomCount(set); i++) {
            expect(Atom.areEqual(AtomSet.atomGetAt(set, i), ret[i])).toBe(true);
        }
    });

    it('element at / index of', () => {
        const control: Atom[] = [];
        const sets = Object.create(null);
        for (let i = 1; i < 10; i++) {
            const set = [];
            for (let j = 1; j < 7; j++) {
                control[control.length] = p(i * i, j * j + 1);
                set[set.length] = j * j + 1;
            }
            sets[i * i] = OrderedSet.ofSortedArray(set);
        }
        const ms = AtomSet.create(sets);
        for (let i = 0; i < control.length; i++) {
            expect(Atom.areEqual(AtomSet.atomGetAt(ms, i), control[i])).toBe(true);
        }

        for (let i = 0; i < control.length; i++) {
            expect(AtomSet.atomIndexOf(ms, control[i])).toBe(i);
        }
    });

    it('packed pairs', () => {
        const set = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        expect(setToPairs(set)).toEqual([p(0, 1), p(0, 2), p(0, 6), p(1, 3)]);
    });

    it('equality', () => {
        const a = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const b = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const c = AtomSet.create([p(1, 3), p(0, 4), p(0, 6), p(0, 2)]);
        const d = AtomSet.create([p(1, 3)]);
        const e = AtomSet.create([p(1, 3)]);
        const f = AtomSet.create([p(3, 3)]);

        expect(AtomSet.areEqual(a, a)).toBe(true);
        expect(AtomSet.areEqual(a, b)).toBe(true);
        expect(AtomSet.areEqual(a, c)).toBe(false);
        expect(AtomSet.areEqual(a, d)).toBe(false);
        expect(AtomSet.areEqual(d, d)).toBe(true);
        expect(AtomSet.areEqual(d, e)).toBe(true);
        expect(AtomSet.areEqual(d, f)).toBe(false);
    });

    it('are intersecting', () => {
        const a = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const b = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const c = AtomSet.create([p(1, 3), p(0, 4), p(0, 6), p(0, 2)]);
        const d = AtomSet.create([p(1, 3)]);
        const e = AtomSet.create([p(1, 3)]);
        const f = AtomSet.create([p(3, 3)]);
        const g = AtomSet.create([p(10, 3), p(8, 1), p(7, 6), p(3, 2)]);

        expect(AtomSet.areIntersecting(a, a)).toBe(true);
        expect(AtomSet.areIntersecting(a, b)).toBe(true);
        expect(AtomSet.areIntersecting(a, c)).toBe(true);
        expect(AtomSet.areIntersecting(a, d)).toBe(true);
        expect(AtomSet.areIntersecting(a, g)).toBe(false);
        expect(AtomSet.areIntersecting(d, d)).toBe(true);
        expect(AtomSet.areIntersecting(d, e)).toBe(true);
        expect(AtomSet.areIntersecting(d, f)).toBe(false);
    });

    it('intersection', () => {
        const a = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const b = AtomSet.create([p(10, 3), p(0, 1), p(0, 6), p(4, 2)]);
        const c = AtomSet.create([p(1, 3)]);
        const d = AtomSet.create([p(2, 3)]);
        expect(AtomSet.intersect(a, a)).toBe(a);
        expect(setToPairs(AtomSet.intersect(a, b))).toEqual([p(0, 1), p(0, 6)]);
        expect(setToPairs(AtomSet.intersect(a, c))).toEqual([p(1, 3)]);
        expect(setToPairs(AtomSet.intersect(c, d))).toEqual([]);
    });

    it('subtract', () => {
        const a = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const a1 = AtomSet.create([p(1, 3), p(0, 1), p(0, 6), p(0, 2)]);
        const b = AtomSet.create([p(10, 3), p(0, 1), p(0, 6), p(4, 2)]);
        const c = AtomSet.create([p(1, 3)]);
        const d = AtomSet.create([p(2, 3)]);
        const e = AtomSet.create([p(0, 2)]);
        expect(setToPairs(AtomSet.subtract(a, a))).toEqual([]);
        expect(setToPairs(AtomSet.subtract(a, a1))).toEqual([]);
        expect(setToPairs(AtomSet.subtract(a, b))).toEqual([p(0, 2), p(1, 3)]);
        expect(setToPairs(AtomSet.subtract(c, d))).toEqual([p(1, 3)]);
        expect(setToPairs(AtomSet.subtract(a, c))).toEqual([p(0, 1), p(0, 2), p(0, 6)]);
        expect(setToPairs(AtomSet.subtract(c, a))).toEqual([]);
        expect(setToPairs(AtomSet.subtract(d, a))).toEqual([p(2, 3)]);
        expect(setToPairs(AtomSet.subtract(a, e))).toEqual([p(0, 1), p(0, 6), p(1, 3)]);
    });

    it('union', () => {
        const a = AtomSet.create([p(1, 3), p(0, 1)]);
        const a1 = AtomSet.create([p(1, 3), p(0, 1)]);
        const b = AtomSet.create([p(10, 3), p(0, 1)]);
        const c = AtomSet.create([p(1, 3)]);
        const d = AtomSet.create([p(2, 3)]);
        expect(AtomSet.unionMany([a])).toBe(a);
        expect(AtomSet.union(a, a)).toBe(a);
        expect(setToPairs(AtomSet.union(a, a))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(AtomSet.union(a, a1))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(AtomSet.union(a, b))).toEqual([p(0, 1), p(1, 3), p(10, 3)]);
        expect(setToPairs(AtomSet.union(c, d))).toEqual([p(1, 3), p(2, 3)]);
        expect(setToPairs(AtomSet.unionMany([a, b, c, d]))).toEqual([p(0, 1), p(1, 3), p(2, 3), p(10, 3)]);
    });
});