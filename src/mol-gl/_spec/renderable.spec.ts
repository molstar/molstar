/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { calculateBoundingSphere } from '../renderable/util';

describe('renderable', () => {
    it('calculateBoundingSphere', () => {
        const position = new Float32Array([
            0, 0, 0,
            1, 0, 0
        ])
        const transform = new Float32Array([
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 0,

            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            1, 0, 0, 0,

            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            2, 0, 0, 0
        ])

        const { boundingSphere } = calculateBoundingSphere(
            position, position.length / 3,
            transform, transform.length / 16
        )

        expect(boundingSphere.radius).toBe(1.5)
        expect(boundingSphere.center).toEqual([1.5, 0.0, 0.0])
    })
})
