/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { EquivalenceClasses } from '../util';

describe('equiv-classes', () => {
    it('integer mod classes', () => {
        const cls = EquivalenceClasses<number, number>(x => x % 2, (a, b) => (a - b) % 2 === 0);
        for (let i = 0; i < 6; i++) cls.add(i, i);

        expect(cls.groups.length).toBe(2);
        expect(cls.groups[0]).toEqual([0, 2, 4]);
        expect(cls.groups[1]).toEqual([1, 3, 5]);
    });
});