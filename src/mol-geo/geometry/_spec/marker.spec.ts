/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { getMarkersAverage } from '../marker-data';

describe('marker-data', () => {
    it('getMarkersAverage', () => {
        expect(getMarkersAverage(new Uint8Array([0, 0, 0, 0]), 3)).toBe(0);
        expect(getMarkersAverage(new Uint8Array([0, 0, 1, 0]), 3)).toBe(1 / 3);
        expect(getMarkersAverage(new Uint8Array([0, 0, 0, 0]), 4)).toBe(0);
        expect(getMarkersAverage(new Uint8Array([0, 0, 1, 0]), 4)).toBe(1 / 4);
    });
});
