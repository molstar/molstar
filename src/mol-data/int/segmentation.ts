/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Interval from './interval';
import OrderedSet from './ordered-set';
import * as Impl from './impl/segmentation';

namespace Segmentation {
    export interface Segment<I extends number = number> { index: I, start: number, end: number }

    export const create: <T extends number = number, I extends number = number>(segs: ArrayLike<T>) => Segmentation<T, I> = Impl.create as any;
    export const ofOffsets: <T extends number = number, I extends number = number>(offsets: ArrayLike<T>, bounds: Interval) => Segmentation<T, I> = Impl.ofOffsets as any;

    export const count: <T extends number = number, I extends number = number>(segs: Segmentation<T, I>) => number = Impl.count as any;
    export const getSegment: <T extends number = number, I extends number = number>(segs: Segmentation<T, I>, value: T) => number = Impl.getSegment as any;
    export const projectValue: <T extends number = number, I extends number = number>(segs: Segmentation<T, I>, set: OrderedSet<T>, value: T) => Interval = Impl.projectValue as any;

    /** Segment iterator that mutates a single segment object to mark all the segments. */
    export const transientSegments: <T extends number = number, I extends number = number>(segs: Segmentation<T, I>, set: OrderedSet<T>, segment?: Segment) => Impl.SegmentIterator<I> = Impl.segments as any;

    export type SegmentIterator<I extends number = number> = Impl.SegmentIterator<I>
}

interface Segmentation<T extends number = number, I extends number = number> {
    '@type': 'segmentation',
    /** All segments are defined by offsets [offsets[i], offsets[i + 1]) for i \in [0, count - 1] */
    readonly offsets: ArrayLike<T>,
    /** Segment index of the i-th element */
    readonly index: ArrayLike<I>,
    readonly count: number
}

export default Segmentation;