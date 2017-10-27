/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/sorted-array'
import Interval from './interval'

namespace SortedArray {
    /** Create interval [min, max] */
    export const create: (xs: ArrayLike<number>) => SortedArray = Impl.ofUnsortedArray as any;
    /** Create interval [min, max) */
    export const ofSortedArray: (xs: ArrayLike<number>) => SortedArray = Impl.ofSortedArray as any;
    export const is: (v: any) => v is Interval = Impl.is as any;

    export const has: (interval: SortedArray, x: number) => boolean = Impl.has as any;
    export const indexOf: (interval: SortedArray, x: number) => number = Impl.indexOf as any;
    export const getAt: (interval: SortedArray, i: number) => number = Impl.getAt as any;

    export const start: (interval: SortedArray) => number = Impl.start as any;
    export const end: (interval: SortedArray) => number = Impl.end as any;
    export const min: (interval: SortedArray) => number = Impl.min as any;
    export const max: (interval: SortedArray) => number = Impl.max as any;
    export const size: (interval: SortedArray) => number = Impl.size as any;
    export const hashCode: (interval: SortedArray) => number = Impl.hashCode as any;

    export const areEqual: (a: SortedArray, b: SortedArray) => boolean = Impl.areEqual as any;

    export const findPredecessorIndex: (interval: SortedArray, x: number) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: (interval: SortedArray, x: number, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: (interval: SortedArray, min: number, max: number) => Interval = Impl.findRange as any;
}

interface SortedArray { '@type': 'int-sorted-array' }

export default SortedArray