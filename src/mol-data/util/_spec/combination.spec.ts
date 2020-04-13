/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { combinations } from '../combination';

describe('Combination', () => {
    it('test123-2', () => {
        const c = combinations([1, 2, 3], 2);
        expect(c).toEqual([[1, 2], [1, 3], [2, 3]]);
    });

    it('test1234-2', () => {
        const c = combinations([1, 2, 3, 4], 2);
        expect(c).toEqual([[1, 2], [1, 3], [2, 3], [1, 4], [2, 4], [3, 4]]);
    });

    it('test1234-1', () => {
        const c = combinations([1, 2, 3, 4], 1);
        expect(c).toEqual([[1], [2], [3], [4]]);
    });

    it('test1234-4', () => {
        const c = combinations([1, 2, 3, 4], 4);
        expect(c).toEqual([[1, 2, 3, 4]]);
    });
});