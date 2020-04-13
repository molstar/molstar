/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import IntTuple from '../tuple';

describe('int pair', () => {
    it('works', () => {
        for (let i = 0; i < 10; i++) {
            for (let j = -10; j < 5; j++) {
                const t = IntTuple.create(i, j);
                expect(IntTuple.fst(t)).toBe(i);
                expect(IntTuple.snd(t)).toBe(j);
            }
        }
    });
});