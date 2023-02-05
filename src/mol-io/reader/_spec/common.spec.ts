/**
 * Copyright (c) 2019-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseFloat as fastParseFloat, parseInt as fastParseInt, getNumberType, NumberTypes } from '../../../mol-io/reader/common/text/number-parser';

describe('common', () => {
    it('number-parser fastParseFloat', () => {
        // ignore suffix numbers in parentheses
        expect(fastParseFloat('11.0829(23)', 0, 11)).toBe(11.0829);
        // scientific with no space between consecutive values
        expect(fastParseFloat('-5.1E-01-6.1E-01', 0, 11)).toBe(-0.51);
        // ignore plus sign
        expect(fastParseFloat('+0.149', 0, 6)).toBe(0.149);
    });

    it('number-parser fastParseInt', () => {
        // ignore suffix numbers in parentheses
        expect(fastParseInt('11(23)', 0, 11)).toBe(11);
        // ignore plus sign
        expect(fastParseFloat('+10149', 0, 6)).toBe(10149);
    });

    it('number-parser getNumberType', () => {
        expect(getNumberType('11')).toBe(NumberTypes.Int);
        expect(getNumberType('5E93')).toBe(NumberTypes.Scientific);
        expect(getNumberType('0.42')).toBe(NumberTypes.Float);
        expect(getNumberType('Foo123')).toBe(NumberTypes.NaN);
        expect(getNumberType('11.0829(23)')).toBe(NumberTypes.NaN);
        expect(getNumberType('1..2')).toBe(NumberTypes.NaN);
        expect(getNumberType('.')).toBe(NumberTypes.NaN);
        expect(getNumberType('-.')).toBe(NumberTypes.NaN);
        expect(getNumberType('e')).toBe(NumberTypes.NaN);
        expect(getNumberType('-e')).toBe(NumberTypes.NaN);
        expect(getNumberType('1e')).toBe(NumberTypes.Scientific);
        expect(getNumberType('-1e')).toBe(NumberTypes.Scientific);
    });
});