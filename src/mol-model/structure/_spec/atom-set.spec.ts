/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int'
import AtomSet from '../structure/atom/set'
import Atom from '../structure/atom'
import AtomGroup from '../structure/atom/group'

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
        const set = AtomSet.ofAtoms([p(10, 11)], AtomSet.Empty);
        expect(setToPairs(set)).toEqual([p(10, 11)]);
        expect(AtomSet.atomHas(set, p(10, 11))).toBe(true);
        expect(AtomSet.atomHas(set, p(11, 11))).toBe(false);
        expect(AtomSet.atomGetAt(set, 0)).toBe(p(10, 11));
        expect(AtomSet.atomCount(set)).toBe(1);
    });

    it('multi', () => {
        const gen = AtomSet.Generator();
        gen.add(1, AtomGroup.createNew(OrderedSet.ofSortedArray([4, 6, 7])));
        gen.add(3, AtomGroup.createNew(OrderedSet.ofRange(0, 1)));
        const set = gen.getSet();
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

    it('template', () => {
        const template = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty)
        const gen = AtomSet.TemplateGenerator(template);
        gen.add(0, AtomGroup.createNew(OrderedSet.ofSortedArray([1, 2, 6])));
        gen.add(1, AtomGroup.createNew(OrderedSet.ofSingleton(3)));
        const set = gen.getSet();

        expect(AtomSet.unitGetById(set, 0)).toBe(AtomSet.unitGetById(template, 0));
        expect(AtomSet.unitGetById(set, 1)).toBe(AtomSet.unitGetById(template, 1));
        expect(set).toBe(template);
    });

    it('template 1', () => {
        const template = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty)
        const gen = AtomSet.TemplateGenerator(template);
        gen.add(0, AtomGroup.createNew(OrderedSet.ofSortedArray([1, 2, 6])));
        gen.add(1, AtomGroup.createNew(OrderedSet.ofSingleton(4)));
        const set = gen.getSet();

        expect(AtomSet.unitGetById(set, 0)).toBe(AtomSet.unitGetById(template, 0));
        expect(AtomSet.unitGetById(set, 1) === AtomSet.unitGetById(template, 1)).toBe(false);
        expect(set === template).toBe(false);
    });

    it('template union', () => {
        const template = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty)

        const p13 = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const p01 = AtomSet.ofAtoms([p(0, 1)], AtomSet.Empty);
        const p02 = AtomSet.ofAtoms([p(0, 2)], AtomSet.Empty);
        const p06 = AtomSet.ofAtoms([p(0, 6)], AtomSet.Empty);

        const u0 = AtomSet.union([p01, p02, p06], template);
        const u1 = AtomSet.union([p01, p02, p06, p13], template);
        expect(AtomSet.unitGetById(u0, 0)).toBe(AtomSet.unitGetById(template, 0));
        expect(AtomSet.unitGetById(u1, 0)).toBe(AtomSet.unitGetById(template, 0));
        expect(AtomSet.unitGetById(u1, 1)).toBe(AtomSet.unitGetById(template, 1));
        expect(u1).toBe(template);
    });

    it('element at / index of', () => {
        const control: Atom[] = [];
        const gen = AtomSet.Generator();
        for (let i = 1; i < 10; i++) {
            const set = [];
            for (let j = 1; j < 7; j++) {
                control[control.length] = p(i * i, j * j + 1);
                set[set.length] = j * j + 1;
            }
            gen.add(i * i, AtomGroup.createNew(OrderedSet.ofSortedArray(set)));
        }
        const ms = gen.getSet();
        for (let i = 0; i < control.length; i++) {
            expect(Atom.areEqual(AtomSet.atomGetAt(ms, i), control[i])).toBe(true);
        }

        for (let i = 0; i < control.length; i++) {
            expect(AtomSet.atomIndexOf(ms, control[i])).toBe(i);
        }
    });

    it('packed pairs', () => {
        const set = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        expect(setToPairs(set)).toEqual([p(0, 1), p(0, 2), p(0, 6), p(1, 3)]);
    });

    it('equality', () => {
        const a = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const b = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const c = AtomSet.ofAtoms([p(1, 3), p(0, 4), p(0, 6), p(0, 2)], AtomSet.Empty);
        const d = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const e = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const f = AtomSet.ofAtoms([p(3, 3)], AtomSet.Empty);

        expect(AtomSet.areEqual(a, a)).toBe(true);
        expect(AtomSet.areEqual(a, b)).toBe(true);
        expect(AtomSet.areEqual(a, c)).toBe(false);
        expect(AtomSet.areEqual(a, d)).toBe(false);
        expect(AtomSet.areEqual(d, d)).toBe(true);
        expect(AtomSet.areEqual(d, e)).toBe(true);
        expect(AtomSet.areEqual(d, f)).toBe(false);
    });

    it('are intersecting', () => {
        const a = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const b = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const c = AtomSet.ofAtoms([p(1, 3), p(0, 4), p(0, 6), p(0, 2)], AtomSet.Empty);
        const d = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const e = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const f = AtomSet.ofAtoms([p(3, 3)], AtomSet.Empty);
        const g = AtomSet.ofAtoms([p(10, 3), p(8, 1), p(7, 6), p(3, 2)], AtomSet.Empty);

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
        const a = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const b = AtomSet.ofAtoms([p(10, 3), p(0, 1), p(0, 6), p(4, 2)], AtomSet.Empty);
        const c = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const d = AtomSet.ofAtoms([p(2, 3)], AtomSet.Empty);
        expect(AtomSet.intersect(a, a)).toBe(a);
        expect(setToPairs(AtomSet.intersect(a, b))).toEqual([p(0, 1), p(0, 6)]);
        expect(setToPairs(AtomSet.intersect(a, c))).toEqual([p(1, 3)]);
        expect(setToPairs(AtomSet.intersect(c, d))).toEqual([]);
    });

    it('subtract', () => {
        const a = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const a1 = AtomSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], AtomSet.Empty);
        const b = AtomSet.ofAtoms([p(10, 3), p(0, 1), p(0, 6), p(4, 2)], AtomSet.Empty);
        const c = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const d = AtomSet.ofAtoms([p(2, 3)], AtomSet.Empty);
        const e = AtomSet.ofAtoms([p(0, 2)], AtomSet.Empty);
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
        const a = AtomSet.ofAtoms([p(1, 3), p(0, 1)], AtomSet.Empty);
        const a1 = AtomSet.ofAtoms([p(1, 3), p(0, 1)], AtomSet.Empty);
        const b = AtomSet.ofAtoms([p(10, 3), p(0, 1)], AtomSet.Empty);
        const c = AtomSet.ofAtoms([p(1, 3)], AtomSet.Empty);
        const d = AtomSet.ofAtoms([p(2, 3)], AtomSet.Empty);
        expect(AtomSet.union([a], AtomSet.Empty)).toBe(a);
        expect(AtomSet.union([a, a], AtomSet.Empty)).toBe(a);
        expect(setToPairs(AtomSet.union([a, a], AtomSet.Empty))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(AtomSet.union([a, a1], AtomSet.Empty))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(AtomSet.union([a, b], AtomSet.Empty))).toEqual([p(0, 1), p(1, 3), p(10, 3)]);
        expect(setToPairs(AtomSet.union([c, d], AtomSet.Empty))).toEqual([p(1, 3), p(2, 3)]);
        expect(setToPairs(AtomSet.union([a, b, c, d], AtomSet.Empty))).toEqual([p(0, 1), p(1, 3), p(2, 3), p(10, 3)]);
    });
});