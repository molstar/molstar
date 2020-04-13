/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat3, Vec3 } from '../3d';

describe('Mat3', () => {
    it('symmetricEigenvalues', () => {
        const m = Mat3.create(
            0.1945, -0.0219, -0.0416,
            -0.0219, 0.1995, -0.0119,
            -0.0416, -0.0119, 0.3673
        );
        const e = Vec3.create(0.377052701425898, 0.21713981522725134, 0.1671074833468507);
        expect(Vec3.equals(e, Mat3.symmetricEigenvalues(Vec3(), m))).toBe(true);
    });

    it('eigenvectors', () => {
        const m = Mat3.create(
            0.1945, -0.0219, -0.0416,
            -0.0219, 0.1995, -0.0119,
            -0.0416, -0.0119, 0.3673
        );
        const e = Vec3.create(0.377052701425898, 0.21713981522725134, 0.1671074833468507);
        const v0 = Vec3.create(-0.2176231019882068, -0.038522620041966125, 0.9752723687391808);
        const v1 = Vec3.create(-0.5905636938047126, 0.8007524989198634, -0.10014968314142503);
        const v2 = Vec3.create(0.7770937582036648, 0.5977553372576602, 0.19701230352667118);

        expect(Vec3.equals(v0, Mat3.eigenvector(Vec3(), m, e[0]))).toBe(true);
        expect(Vec3.equals(v1, Mat3.eigenvector(Vec3(), m, e[1]))).toBe(true);
        expect(Vec3.equals(v2, Mat3.eigenvector(Vec3(), m, e[2]))).toBe(true);
    });
});