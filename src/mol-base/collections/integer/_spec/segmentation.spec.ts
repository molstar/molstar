/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import OrderedSet from '../ordered-set'
import Interval from '../interval'
import Segmentation from '../segmentation'

describe('segments', () => {
    const data = OrderedSet.ofSortedArray([4, 9, 10, 11, 14, 15, 16]);
    const segs = Segmentation.create([0, 4, 10, 12, 13, 15, 25]);

    it('size', () => expect(Segmentation.count(segs)).toBe(6));

    it('project', () => {
        const p = Segmentation.projectValue(segs, data, 4);
        expect(p).toBe(Interval.ofBounds(0, 2))
    });

    it('ofOffsetts', () => {
        const p = Segmentation.ofOffsets([10, 12], Interval.ofBounds(10, 14));
        expect(p.segments).toEqual(new Int32Array([0, 2, 4]))
    });

    it('map', () => {
        const segs = Segmentation.create([0, 1, 2]);
        expect(segs.segmentMap).toEqual(new Int32Array([0, 1]));
        expect(Segmentation.getSegment(segs, 0)).toBe(0);
        expect(Segmentation.getSegment(segs, 1)).toBe(1);
    });

    it('iteration', () => {
        const it = Segmentation.transientSegments(segs, data);

        const t = Object.create(null);
        while (it.hasNext) {
            const s = it.move();
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
