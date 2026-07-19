/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Spacegroup } from '../../../../mol-math/geometry';
import { RuntimeContext } from '../../../../mol-task';
import { computeElectronDensityFromReflections } from '../reflections';

describe('computeElectronDensityFromReflections', () => {
    it('sizes the reciprocal grid from symmetry-expanded reflections', async () => {
        const swapHK = [
            0, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
        ];
        const sg = { operators: [swapHK] } as unknown as Spacegroup;

        const { N0, N1, N2 } = await computeElectronDensityFromReflections(
            Int16Array.of(1),
            Int16Array.of(9),
            Int16Array.of(0),
            Float32Array.of(1),
            Float32Array.of(0),
            sg,
            RuntimeContext.Synchronous,
        );

        expect(N0).toBe(32);
        expect(N1).toBe(4);
        expect(N2).toBe(2);
    });
});