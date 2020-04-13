/**
 * Copyright (c) 2017-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/sorted-array';
import Interval from './interval';

namespace SortedArray {
    export const Empty: SortedArray = Impl.Empty as any;
    export const ofUnsortedArray: <T extends number = number>(xs: ArrayLike<number>) => SortedArray<T> = Impl.ofUnsortedArray as any;
    export const ofSingleton: <T extends number = number>(v: number) => SortedArray<T> = Impl.ofSingleton as any;
    export const ofSortedArray: <T extends number = number>(xs: ArrayLike<number>) => SortedArray<T> = Impl.ofSortedArray as any;
    /** create sorted array [min, max] (it DOES contain the max value) */
    export const ofRange: <T extends number = number>(min: T, max: T) => SortedArray<T> = Impl.ofRange as any;
    /** create sorted array [min, max) (it does NOT contain the max value) */
    export const ofBounds: <T extends number = number>(min: T, max: T) => SortedArray<T> = (min, max) => Impl.ofRange(min, max - 1) as any;
    export const is: <T extends number = number>(v: any) => v is SortedArray<T> = Impl.is as any;

    export const has: <T extends number = number>(array: SortedArray<T>, x: T) => boolean = Impl.has as any;
    /** Returns the index of `x` in `set` or -1 if not found. */
    export const indexOf: <T extends number = number>(array: SortedArray<T>, x: T) => number = Impl.indexOf as any;
    export const indexOfInInterval: <T extends number = number>(array: SortedArray<T>, x: number, bounds: Interval) => number = Impl.indexOfInInterval as any;
    export const indexOfInRange: <T extends number = number>(array: SortedArray<T>, x: number, start: number, end: number) => number = Impl.indexOfInRange as any;

    /** Returns `array[0]` */
    export const start: <T extends number = number>(array: SortedArray<T>) => T = Impl.start as any;
    /** Returns `array[array.length - 1] + 1` */
    export const end: <T extends number = number>(array: SortedArray<T>) => T = Impl.end as any;
    export const min: <T extends number = number>(array: SortedArray<T>) => T = Impl.min as any;
    export const max: <T extends number = number>(array: SortedArray<T>) => T = Impl.max as any;
    export const size: <T extends number = number>(array: SortedArray<T>) => number = Impl.size as any;
    export const hashCode: <T extends number = number>(array: SortedArray<T>) => number = Impl.hashCode as any;
    export const toString: <T extends number = number>(array: SortedArray<T>) => string = Impl.toString as any;

    export const areEqual: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => boolean = Impl.areEqual as any;
    export const areIntersecting: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => boolean = Impl.areIntersecting as any;
    export const isSubset: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => boolean = Impl.isSubset as any;

    export const union: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => SortedArray<T> = Impl.union as any;
    export const intersect: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => SortedArray<T> = Impl.intersect as any;
    export const subtract: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => SortedArray<T> = Impl.subtract as any;

    export const findPredecessorIndex: <T extends number = number>(array: SortedArray<T>, x: T) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: <T extends number = number>(array: SortedArray<T>, x: T, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: <T extends number = number>(array: SortedArray<T>, min: T, max: T) => Interval = Impl.findRange as any;
    export const intersectionSize: <T extends number = number>(a: SortedArray<T>, b: SortedArray<T>) => number = Impl.intersectionSize as any;

    export const deduplicate: <T extends number = number>(array: SortedArray<T>) => SortedArray<T> = Impl.deduplicate as any;
    /** Returns indices of xs in the array. E.g. indicesOf([10, 11, 12], [10, 12]) ==> [0, 2] */
    export const indicesOf: <T extends number = number, I extends number = number>(array: SortedArray<T>, xs: SortedArray<T>) => SortedArray<I> = Impl.indicesOf as any;
}

interface SortedArray<T extends number = number> extends ArrayLike<T> { '@type': 'int-sorted-array' }

export default SortedArray;