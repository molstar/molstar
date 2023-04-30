/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { NumberArray } from '../../../mol-util/type-helpers';
import { PrincipalAxes } from '../matrix/principal-axes';

describe('PrincipalAxes', () => {
    it('same-cartesian-plane', () => {
        const positions: NumberArray = [ // same y coordinate
            0.1945, -0.0219, -0.0416,
            -0.0219, -0.0219, -0.0119,
        ];
        const { origin } = PrincipalAxes.ofPositions(positions).boxAxes;
        expect(origin[0] !== Infinity && origin[1] !== Infinity && origin[2] !== Infinity).toBe(true);
    });
});
