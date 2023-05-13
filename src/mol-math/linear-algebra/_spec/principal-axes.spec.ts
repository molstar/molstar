/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { NumberArray } from '../../../mol-util/type-helpers';
import { Vec3 } from '../3d/vec3';
import { PrincipalAxes } from '../matrix/principal-axes';

describe('PrincipalAxes', () => {
    it('same-cartesian-plane', () => {
        const positions: NumberArray = [ // same y coordinate
            0.1945, -0.0219, -0.0416,
            -0.0219, -0.0219, -0.0119,
        ];
        const pa = PrincipalAxes.ofPositions(positions);
        expect(Vec3.isFinite(pa.boxAxes.origin)).toBe(true);
        expect(Vec3.equals(pa.boxAxes.origin, pa.momentsAxes.origin)).toBe(true);
    });

    it('same-point', () => {
        const positions: NumberArray = [ // same coordinates
            0.1945, -0.0219, -0.0416,
            0.1945, -0.0219, -0.0416,
        ];
        const pa = PrincipalAxes.ofPositions(positions);
        expect(Vec3.isFinite(pa.boxAxes.origin)).toBe(true);
        expect(Vec3.equals(pa.boxAxes.origin, pa.momentsAxes.origin)).toBe(true);
    });
});
