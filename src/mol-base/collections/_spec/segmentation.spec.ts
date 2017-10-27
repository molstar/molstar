/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from '../ordered-set'
import Interval from '../interval'
import Segmentation from '../segmentation'
import SortedArray from '../sorted-array'

describe('segments', () => {
    const data = OrderedSet.ofSortedArray([4, 9, 10, 11, 14, 15, 16]);
    const segs = Segmentation.create(SortedArray.ofSortedArray([0, 4, 10, 12, 13, 15, 25]), [])

    it('project', () => {
        const p = Segmentation.projectValue(segs, data, 4);
        expect(p).toBe(Interval.ofBounds(0, 2))
    });

    it('iteration', () => {
        const it = Segmentation.segments(segs, data);

        const t = Object.create(null);
        for (let s = it.move(); !it.done; s = it.move()) {
            for (let j = s.start; j < s.end; j++) {
                const x = t[s.index];
                const v = OrderedSet.getAt(data, j);
                if (!x) t[s.index] = [v];
                else x[x.length] = v;
            }
        }

        expect(t).toEqual({ 1: [4, 9], 2: [10, 11], 4: [14], 5: [15, 16] });
    });
});
