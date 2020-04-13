/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseFloat as fastParseFloat, parseInt as fastParseInt, getNumberType, NumberType } from '../../../mol-io/reader/common/text/number-parser';

describe('common', () => {
    it('number-parser fastParseFloat', () => {
        expect(fastParseFloat('11.0829(23)', 0, 11)).toBe(11.0829);
    });

    it('number-parser fastParseInt', () => {
        expect(fastParseInt('11(23)', 0, 11)).toBe(11);
    });

    it('number-parser getNumberType', () => {
        expect(getNumberType('11')).toBe(NumberType.Int);
        expect(getNumberType('5E93')).toBe(NumberType.Scientific);
        expect(getNumberType('0.42')).toBe(NumberType.Float);
        expect(getNumberType('Foo123')).toBe(NumberType.NaN);
        expect(getNumberType('11.0829(23)')).toBe(NumberType.NaN);
        expect(getNumberType('1..2')).toBe(NumberType.NaN);
        expect(getNumberType('.')).toBe(NumberType.NaN);
        expect(getNumberType('-.')).toBe(NumberType.NaN);
        expect(getNumberType('e')).toBe(NumberType.NaN);
        expect(getNumberType('-e')).toBe(NumberType.NaN);
        expect(getNumberType('1e')).toBe(NumberType.Scientific);
        expect(getNumberType('-1e')).toBe(NumberType.Scientific);
    });
});