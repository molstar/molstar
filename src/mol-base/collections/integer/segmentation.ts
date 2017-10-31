/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../iterator'
import Interval from './interval'
import OrderedSet from './ordered-set'
import * as Impl from './impl/segmentation'

namespace Segmentation {
    export interface Segment { index: number, start: number, end: number }

    export const create: (segs: ArrayLike<number>) => Segmentation = Impl.create as any;
    export const ofOffsets: (offsets: ArrayLike<number>, bounds: Interval) => Segmentation = Impl.ofOffsets as any;

    export const count: (segs: Segmentation) => number = Impl.count as any;
    export const getSegment: (segs: Segmentation, value: number) => number = Impl.getSegment as any;
    export const projectValue: (segs: Segmentation, set: OrderedSet, value: number) => Interval = Impl.projectValue as any;

    // Segment iterator that mutates a single segment object to mark all the segments.
    export const transientSegments: (segs: Segmentation, set: OrderedSet, segment?: Segment) => Iterator<Segment> = Impl.segments as any;
}

interface Segmentation {
    '@type': 'segmentation',
    readonly segments: ArrayLike<number>,
    readonly segmentMap: ArrayLike<number>,
    readonly count: number
}

export default Segmentation