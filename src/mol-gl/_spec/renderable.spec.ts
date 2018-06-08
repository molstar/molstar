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

        const bs = calculateBoundingSphere({
            position,
            positionCount: position.length / 3,
            transform,
            transformCount: transform.length / 16
        })

        expect(bs.radius).toBe(1.5)
        expect(bs.center).toEqual([1.5, 0.0, 0.0])
    })
})
