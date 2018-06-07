/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/sorted-array'
import Interval from './interval'

namespace SortedArray {
    export const Empty: SortedArray = Impl.Empty as any;
    export const ofUnsortedArray: (xs: ArrayLike<number>) => SortedArray = Impl.ofUnsortedArray as any;
    export const ofSingleton: (v: number) => SortedArray = Impl.ofSingleton as any;
    export const ofSortedArray: (xs: ArrayLike<number>) => SortedArray = Impl.ofSortedArray as any;
    // create sorted array [min, max] (it DOES contain the max value)
    export const ofRange: (min: number, max: number) => SortedArray = Impl.ofRange as any;
    // create sorted array [min, max) (it DOES not contain the max value)
    export const ofBounds: (min: number, max: number) => SortedArray = (min, max) => Impl.ofRange(min, max - 1) as any;
    export const is: (v: any) => v is Interval = Impl.is as any;

    export const has: (array: SortedArray, x: number) => boolean = Impl.has as any;
    export const indexOf: (array: SortedArray, x: number) => number = Impl.indexOf as any;
    export const indexOfInInterval: (array: SortedArray, x: number, bounds: Interval) => number = Impl.indexOfInInterval as any;

    export const start: (array: SortedArray) => number = Impl.start as any;
    export const end: (array: SortedArray) => number = Impl.end as any;
    export const min: (array: SortedArray) => number = Impl.min as any;
    export const max: (array: SortedArray) => number = Impl.max as any;
    export const size: (array: SortedArray) => number = Impl.size as any;
    export const hashCode: (array: SortedArray) => number = Impl.hashCode as any;

    export const areEqual: (a: SortedArray, b: SortedArray) => boolean = Impl.areEqual as any;
    export const areIntersecting: (a: SortedArray, b: SortedArray) => boolean = Impl.areIntersecting as any;
    export const isSubset: (a: SortedArray, b: SortedArray) => boolean = Impl.isSubset as any;

    export const union: (a: SortedArray, b: SortedArray) => SortedArray = Impl.union as any;
    export const intersect: (a: SortedArray, b: SortedArray) => SortedArray = Impl.intersect as any;
    export const subtract: (a: SortedArray, b: SortedArray) => SortedArray = Impl.subtract as any;

    export const findPredecessorIndex: (array: SortedArray, x: number) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: (array: SortedArray, x: number, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: (array: SortedArray, min: number, max: number) => Interval = Impl.findRange as any;

    export const deduplicate: (array: SortedArray) => SortedArray = Impl.deduplicate as any;
    /** Returns indices of xs in the array. E.g. indicesOf([10, 11, 12], [10, 12]) ==> [0, 2] */
    export const indicesOf: (array: SortedArray, xs: SortedArray) => SortedArray = Impl.indicesOf as any;
}

interface SortedArray extends ArrayLike<number> { '@type': 'int-sorted-array' }

export default SortedArray