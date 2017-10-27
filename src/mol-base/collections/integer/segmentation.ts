/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import Iterator from '../iterator'
import Interval from './interval'
import OrderedSet from './ordered-set'
import SortedArray from './sorted-array'
import * as Impl from './impl/segmentation'

namespace Segmentation {
    export interface Segment { index: number, start: number, end: number }

    export const create: (segs: SortedArray, segIndex: ArrayLike<number>) => Segmentation = Impl.create as any;

    export const getSegment: (segs: Segmentation, value: number) => number = Impl.getSegment as any;
    export const projectValue: (segs: Segmentation, set: OrderedSet, value: number) => Interval = Impl.projectValue as any;

    export const segments: (segs: Segmentation, set: OrderedSet, range?: Interval) => Iterator<Segment> = Impl.segments as any;
}

interface Segmentation { '@type': 'segmentation' }

export default Segmentation