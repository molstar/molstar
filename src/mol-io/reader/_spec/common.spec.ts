/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { parseFloat as fastParseFloat, parseInt as fastParseInt } from '../../../mol-io/reader/common/text/number-parser';

describe('common', () => {
    it('number-parser fastParseFloat', () => {
        expect(fastParseFloat('11.0829(23)', 0, 11)).toBe(11.0829)
    });

    it('number-parser fastParseInt', () => {
        expect(fastParseInt('11(23)', 0, 11)).toBe(11)
    });
});