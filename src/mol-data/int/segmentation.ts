/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Interval from './interval'
import OrderedSet from './ordered-set'
import * as Impl from './impl/segmentation'

namespace Segmentation {
    export interface Segment<T extends number = number> { index: number, start: number, end: number }

    export const create: <T extends number = number>(segs: ArrayLike<T>) => Segmentation<T> = Impl.create as any;
    export const ofOffsets: <T extends number = number>(offsets: ArrayLike<T>, bounds: Interval) => Segmentation<T> = Impl.ofOffsets as any;

    export const count: <T extends number = number>(segs: Segmentation<T>) => number = Impl.count as any;
    export const getSegment: <T extends number = number>(segs: Segmentation<T>, value: T) => number = Impl.getSegment as any;
    export const projectValue: <T extends number = number>(segs: Segmentation<T>, set: OrderedSet<T>, value: T) => Interval = Impl.projectValue as any;

    // Segment iterator that mutates a single segment object to mark all the segments.
    export const transientSegments: <T extends number = number>(segs: Segmentation<T>, set: OrderedSet<T>, segment?: Segment<T>) => Impl.SegmentIterator = Impl.segments as any;
}

interface Segmentation<T extends number = number> {
    '@type': 'segmentation',
    readonly segments: ArrayLike<T>,
    readonly segmentMap: ArrayLike<number>,
    readonly count: number
}

export default Segmentation