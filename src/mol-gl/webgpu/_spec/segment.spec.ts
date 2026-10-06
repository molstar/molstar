/** Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info. */

import { Tensor } from '../../../mol-math/linear-algebra';
import { createSegmentSampler } from '../../../mol-repr/volume/util';

describe('segment boundary sampling', () => {
    it('zero-pads every grid face without aliasing neighboring rows in any axis order', () => {
        for (const order of [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]) {
            const space = Tensor.Space([3, 4, 5], order, Float32Array), data = space.create(); data.fill(7);
            const sample = createSegmentSampler(Tensor.create(space, data), [7]);
            expect(sample(1, 1, 1)).toBe(255);
            expect(sample(0, 1, 1)).toBe(223);
            expect(sample(0, 0, 0)).toBe(159);
            for (let a = 0; a < 3; a++) for (const side of [-1, 1]) {
                const c = [1, 1, 1]; c[a] = side === -1 ? -1 : space.dimensions[a];
                expect(sample(...c as [number, number, number])).toBe(32);
                c[a] += side;
                expect(sample(...c as [number, number, number])).toBe(0);
            }
            expect(createSegmentSampler(Tensor.create(space, data), [8])(1, 1, 1)).toBe(0);
        }
    });
});
