/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Ccp4Header } from '../../../mol-io/reader/ccp4/schema';
import { Cell } from '../../../mol-math/geometry/spacegroup/cell';
import { Vec3 } from '../../../mol-math/linear-algebra';
import { degToRad } from '../../../mol-math/misc';
import { getCcp4Angles } from '../ccp4';

function headerWithAngles(alpha: number, beta: number, gamma: number) {
    return { alpha, beta, gamma } as unknown as Ccp4Header;
}

describe('getCcp4Angles', () => {
    let warn: jest.SpyInstance;

    beforeEach(() => {
        warn = jest.spyOn(console, 'warn').mockImplementation(() => {});
    });

    afterEach(() => {
        warn.mockRestore();
    });

    it('passes valid angles through as radians', () => {
        const angles = getCcp4Angles(headerWithAngles(90, 90, 120));
        expect(angles[0]).toBeCloseTo(degToRad(90));
        expect(angles[1]).toBeCloseTo(degToRad(90));
        expect(angles[2]).toBeCloseTo(degToRad(120));
        expect(warn).not.toHaveBeenCalled();
    });

    it('assumes a right angle when the angles are unset', () => {
        // MRC volumes written by IMOD leave `cellb` at zero
        const angles = getCcp4Angles(headerWithAngles(0, 0, 0));
        for (let i = 0; i < 3; ++i) expect(angles[i]).toBeCloseTo(degToRad(90));
        expect(warn).toHaveBeenCalledTimes(3);
    });

    it('assumes a right angle for out of range angles', () => {
        const angles = getCcp4Angles(headerWithAngles(180, -5, NaN));
        for (let i = 0; i < 3; ++i) expect(angles[i]).toBeCloseTo(degToRad(90));
    });

    it('gives a usable cell for unset angles', () => {
        // zero angles would make the z2 basis term 0/0, giving a NaN transform
        const size = Vec3.create(6806.36, 9580.8, 3333.32);
        const cell = Cell.create(size, getCcp4Angles(headerWithAngles(0, 0, 0)));
        for (let i = 0; i < 16; ++i) {
            expect(Number.isFinite(cell.fromFractional[i])).toBe(true);
            expect(Number.isFinite(cell.toFractional[i])).toBe(true);
        }
        expect(cell.volume).toBeCloseTo(size[0] * size[1] * size[2]);
    });
});
