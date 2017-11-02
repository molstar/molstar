/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Base from './impl/ordered-set'
import Interval from './interval'

namespace OrderedSet {
    export const Empty: OrderedSet = Base.Empty as any;
    export const ofSingleton: (value: number) => OrderedSet = Base.ofSingleton as any;
    export const ofRange: (min: number, max: number) => OrderedSet = Base.ofRange as any;
    export const ofBounds: (min: number, max: number) => OrderedSet = Base.ofBounds as any;
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export const ofSortedArray: (xs: ArrayLike<number>) => OrderedSet = Base.ofSortedArray as any;

    export const has: (set: OrderedSet, x: number) => boolean = Base.has as any;
    export const indexOf: (set: OrderedSet, x: number) => number = Base.indexOf as any;
    export const getAt: (set: OrderedSet, i: number) => number = Base.getAt as any;

    export const min: (set: OrderedSet) => number = Base.min as any;
    export const max: (set: OrderedSet) => number = Base.max as any;
    export const size: (set: OrderedSet) => number = Base.size as any;
    export const hashCode: (set: OrderedSet) => number = Base.hashCode as any;

    export const areEqual: (a: OrderedSet, b: OrderedSet) => boolean = Base.areEqual as any;
    export const areIntersecting: (a: OrderedSet, b: OrderedSet) => boolean = Base.areIntersecting as any;
    export const isSubset: (a: OrderedSet, b: OrderedSet) => boolean = Base.isSubset as any;

    export const union: (a: OrderedSet, b: OrderedSet) => OrderedSet = Base.union as any;
    export const intersect: (a: OrderedSet, b: OrderedSet) => OrderedSet = Base.intersect as any;
    export const subtract: (a: OrderedSet, b: OrderedSet) => OrderedSet = Base.subtract as any;

    export const findPredecessorIndex: (set: OrderedSet, x: number) => number = Base.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: (set: OrderedSet, x: number, range: Interval) => number = Base.findPredecessorIndexInInterval as any;
    export const findRange: (set: OrderedSet, min: number, max: number) => Interval = Base.findRange as any;
}

interface OrderedSet { '@type': 'int-interval' | 'int-sorted-array' }

export default OrderedSet