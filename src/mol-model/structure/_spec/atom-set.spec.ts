/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { OrderedSet } from 'mol-data/int'
import ElementSet from '../structure/element/set'
import Element from '../structure/element'
import ElementGroup from '../structure/element/group'

describe('atom set', () => {
    const p = (i: number, j: number) => Element.create(i, j);

    function setToPairs(set: ElementSet): ArrayLike<Element> {
        const ret: Element[] = [];
        const it = ElementSet.elements(set);
        while (it.hasNext) {
            ret[ret.length] = it.move();
        }
        return ret;
    }

    it('singleton pair', () => {
        const set = ElementSet.ofAtoms([p(10, 11)], ElementSet.Empty);
        expect(setToPairs(set)).toEqual([p(10, 11)]);
        expect(ElementSet.elementHas(set, p(10, 11))).toBe(true);
        expect(ElementSet.elementHas(set, p(11, 11))).toBe(false);
        expect(ElementSet.elementGetAt(set, 0)).toBe(p(10, 11));
        expect(ElementSet.elementCount(set)).toBe(1);
    });

    it('singleton atom', () => {
        const set = ElementSet.singleton(p(10, 11), ElementSet.Empty);
        expect(setToPairs(set)).toEqual([p(10, 11)]);
        expect(ElementSet.elementHas(set, p(10, 11))).toBe(true);
        expect(ElementSet.elementHas(set, p(11, 11))).toBe(false);
        expect(ElementSet.elementGetAt(set, 0)).toBe(p(10, 11));
        expect(ElementSet.elementCount(set)).toBe(1);
    });

    it('multi', () => {
        const gen = ElementSet.Generator();
        gen.add(1, ElementGroup.createNew(OrderedSet.ofSortedArray([4, 6, 7])));
        gen.add(3, ElementGroup.createNew(OrderedSet.ofRange(0, 1)));
        const set = gen.getSet();
        const ret = [p(1, 4), p(1, 6), p(1, 7), p(3, 0), p(3, 1)];
        expect(ElementSet.elementCount(set)).toBe(ret.length);
        expect(setToPairs(set)).toEqual([p(1, 4), p(1, 6), p(1, 7), p(3, 0), p(3, 1)]);
        expect(ElementSet.elementHas(set, p(10, 11))).toBe(false);
        expect(ElementSet.elementHas(set, p(3, 0))).toBe(true);
        expect(ElementSet.elementHas(set, p(1, 7))).toBe(true);
        for (let i = 0; i < ElementSet.elementCount(set); i++) {
            expect(Element.areEqual(ElementSet.elementGetAt(set, i), ret[i])).toBe(true);
        }
    });

    it('template', () => {
        const template = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty)
        const gen = ElementSet.TemplateGenerator(template);
        gen.add(0, OrderedSet.ofSortedArray([1, 2, 6]));
        gen.add(1, OrderedSet.ofSingleton(3));
        const set = gen.getSet();

        expect(ElementSet.unitGetById(set, 0)).toBe(ElementSet.unitGetById(template, 0));
        expect(ElementSet.unitGetById(set, 1)).toBe(ElementSet.unitGetById(template, 1));
        expect(set).toBe(template);
    });

    it('template 1', () => {
        const template = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty)
        const gen = ElementSet.TemplateGenerator(template);
        gen.add(0, OrderedSet.ofSortedArray([1, 2, 6]));
        gen.add(1, OrderedSet.ofSingleton(4));
        const set = gen.getSet();

        expect(ElementSet.unitGetById(set, 0)).toBe(ElementSet.unitGetById(template, 0));
        expect(ElementSet.unitGetById(set, 1) === ElementSet.unitGetById(template, 1)).toBe(false);
        expect(set === template).toBe(false);
    });

    it('template union', () => {
        const template = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty)

        const p13 = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const p01 = ElementSet.ofAtoms([p(0, 1)], ElementSet.Empty);
        const p02 = ElementSet.ofAtoms([p(0, 2)], ElementSet.Empty);
        const p06 = ElementSet.ofAtoms([p(0, 6)], ElementSet.Empty);

        const u0 = ElementSet.union([p01, p02, p06], template);
        const u1 = ElementSet.union([p01, p02, p06, p13], template);
        expect(ElementSet.unitGetById(u0, 0)).toBe(ElementSet.unitGetById(template, 0));
        expect(ElementSet.unitGetById(u1, 0)).toBe(ElementSet.unitGetById(template, 0));
        expect(ElementSet.unitGetById(u1, 1)).toBe(ElementSet.unitGetById(template, 1));
        expect(u1).toBe(template);
    });

    it('element at / index of', () => {
        const control: Element[] = [];
        const gen = ElementSet.Generator();
        for (let i = 1; i < 10; i++) {
            const set = [];
            for (let j = 1; j < 7; j++) {
                control[control.length] = p(i * i, j * j + 1);
                set[set.length] = j * j + 1;
            }
            gen.add(i * i, ElementGroup.createNew(OrderedSet.ofSortedArray(set)));
        }
        const ms = gen.getSet();
        for (let i = 0; i < control.length; i++) {
            expect(Element.areEqual(ElementSet.elementGetAt(ms, i), control[i])).toBe(true);
        }

        for (let i = 0; i < control.length; i++) {
            expect(ElementSet.elementIndexOf(ms, control[i])).toBe(i);
        }
    });

    it('packed pairs', () => {
        const set = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        expect(setToPairs(set)).toEqual([p(0, 1), p(0, 2), p(0, 6), p(1, 3)]);
    });

    it('equality', () => {
        const a = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const b = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const c = ElementSet.ofAtoms([p(1, 3), p(0, 4), p(0, 6), p(0, 2)], ElementSet.Empty);
        const d = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const e = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const f = ElementSet.ofAtoms([p(3, 3)], ElementSet.Empty);

        expect(ElementSet.areEqual(a, a)).toBe(true);
        expect(ElementSet.areEqual(a, b)).toBe(true);
        expect(ElementSet.areEqual(a, c)).toBe(false);
        expect(ElementSet.areEqual(a, d)).toBe(false);
        expect(ElementSet.areEqual(d, d)).toBe(true);
        expect(ElementSet.areEqual(d, e)).toBe(true);
        expect(ElementSet.areEqual(d, f)).toBe(false);
    });

    it('are intersecting', () => {
        const a = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const b = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const c = ElementSet.ofAtoms([p(1, 3), p(0, 4), p(0, 6), p(0, 2)], ElementSet.Empty);
        const d = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const e = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const f = ElementSet.ofAtoms([p(3, 3)], ElementSet.Empty);
        const g = ElementSet.ofAtoms([p(10, 3), p(8, 1), p(7, 6), p(3, 2)], ElementSet.Empty);

        expect(ElementSet.areIntersecting(a, a)).toBe(true);
        expect(ElementSet.areIntersecting(a, b)).toBe(true);
        expect(ElementSet.areIntersecting(a, c)).toBe(true);
        expect(ElementSet.areIntersecting(a, d)).toBe(true);
        expect(ElementSet.areIntersecting(a, g)).toBe(false);
        expect(ElementSet.areIntersecting(d, d)).toBe(true);
        expect(ElementSet.areIntersecting(d, e)).toBe(true);
        expect(ElementSet.areIntersecting(d, f)).toBe(false);
    });

    it('intersection', () => {
        const a = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const b = ElementSet.ofAtoms([p(10, 3), p(0, 1), p(0, 6), p(4, 2)], ElementSet.Empty);
        const c = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const d = ElementSet.ofAtoms([p(2, 3)], ElementSet.Empty);
        expect(ElementSet.intersect(a, a)).toBe(a);
        expect(setToPairs(ElementSet.intersect(a, b))).toEqual([p(0, 1), p(0, 6)]);
        expect(setToPairs(ElementSet.intersect(a, c))).toEqual([p(1, 3)]);
        expect(setToPairs(ElementSet.intersect(c, d))).toEqual([]);
    });

    it('subtract', () => {
        const a = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const a1 = ElementSet.ofAtoms([p(1, 3), p(0, 1), p(0, 6), p(0, 2)], ElementSet.Empty);
        const b = ElementSet.ofAtoms([p(10, 3), p(0, 1), p(0, 6), p(4, 2)], ElementSet.Empty);
        const c = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const d = ElementSet.ofAtoms([p(2, 3)], ElementSet.Empty);
        const e = ElementSet.ofAtoms([p(0, 2)], ElementSet.Empty);
        expect(setToPairs(ElementSet.subtract(a, a))).toEqual([]);
        expect(setToPairs(ElementSet.subtract(a, a1))).toEqual([]);
        expect(setToPairs(ElementSet.subtract(a, b))).toEqual([p(0, 2), p(1, 3)]);
        expect(setToPairs(ElementSet.subtract(c, d))).toEqual([p(1, 3)]);
        expect(setToPairs(ElementSet.subtract(a, c))).toEqual([p(0, 1), p(0, 2), p(0, 6)]);
        expect(setToPairs(ElementSet.subtract(c, a))).toEqual([]);
        expect(setToPairs(ElementSet.subtract(d, a))).toEqual([p(2, 3)]);
        expect(setToPairs(ElementSet.subtract(a, e))).toEqual([p(0, 1), p(0, 6), p(1, 3)]);
    });

    it('union', () => {
        const a = ElementSet.ofAtoms([p(1, 3), p(0, 1)], ElementSet.Empty);
        const a1 = ElementSet.ofAtoms([p(1, 3), p(0, 1)], ElementSet.Empty);
        const b = ElementSet.ofAtoms([p(10, 3), p(0, 1)], ElementSet.Empty);
        const c = ElementSet.ofAtoms([p(1, 3)], ElementSet.Empty);
        const d = ElementSet.ofAtoms([p(2, 3)], ElementSet.Empty);
        expect(ElementSet.union([a], ElementSet.Empty)).toBe(a);
        expect(ElementSet.union([a, a], ElementSet.Empty)).toBe(a);
        expect(setToPairs(ElementSet.union([a, a], ElementSet.Empty))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(ElementSet.union([a, a1], ElementSet.Empty))).toEqual([p(0, 1), p(1, 3)]);
        expect(setToPairs(ElementSet.union([a, b], ElementSet.Empty))).toEqual([p(0, 1), p(1, 3), p(10, 3)]);
        expect(setToPairs(ElementSet.union([c, d], ElementSet.Empty))).toEqual([p(1, 3), p(2, 3)]);
        expect(setToPairs(ElementSet.union([a, b, c, d], ElementSet.Empty))).toEqual([p(0, 1), p(1, 3), p(2, 3), p(10, 3)]);
    });
});