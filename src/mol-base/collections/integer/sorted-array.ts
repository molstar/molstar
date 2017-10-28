/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/sorted-array'
import Interval from './interval'

namespace SortedArray {
    export const ofUnsortedArray: (xs: ArrayLike<number>) => SortedArray = Impl.ofUnsortedArray as any;
    export const ofSortedArray: (xs: ArrayLike<number>) => SortedArray = Impl.ofSortedArray as any;
    export const is: (v: any) => v is Interval = Impl.is as any;

    export const values: (array: SortedArray) => ArrayLike<number> = (xs) => xs as any;

    export const has: (array: SortedArray, x: number) => boolean = Impl.has as any;
    export const indexOf: (array: SortedArray, x: number) => number = Impl.indexOf as any;
    export const indexOfInterval: (array: SortedArray, x: number, bounds: Interval) => number = Impl.indexOfInterval as any;
    export const getAt: (array: SortedArray, i: number) => number = Impl.getAt as any;

    export const start: (array: SortedArray) => number = Impl.start as any;
    export const end: (array: SortedArray) => number = Impl.end as any;
    export const min: (array: SortedArray) => number = Impl.min as any;
    export const max: (array: SortedArray) => number = Impl.max as any;
    export const size: (array: SortedArray) => number = Impl.size as any;
    export const hashCode: (array: SortedArray) => number = Impl.hashCode as any;

    export const areEqual: (a: SortedArray, b: SortedArray) => boolean = Impl.areEqual as any;

    export const findPredecessorIndex: (array: SortedArray, x: number) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: (array: SortedArray, x: number, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: (array: SortedArray, min: number, max: number) => Interval = Impl.findRange as any;
}

interface SortedArray { '@type': 'int-sorted-array' }

export default SortedArray