/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Base from './ordered-set/base'
import SegmentIterator from './ordered-set/segment-iterator'

namespace OrderedSet {
    export interface IndexRange { start: number, end: number }
    export function IndexRange(start?: number, end?: number): IndexRange { return { start: start || 0, end: end || 0 }; }

    export const Empty: OrderedSet = Base.Empty as any;
    export const ofSingleton: (value: number) => OrderedSet = Base.ofSingleton as any;
    export const ofRange: (min: number, max: number) => OrderedSet = Base.ofRange as any;
    /** It is the responsibility of the caller to ensure the array is sorted and contains unique values. */
    export const ofSortedArray: (xs: ArrayLike<number>) => OrderedSet = Base.ofSortedArray as any;

    export const has: (set: OrderedSet, x: number) => boolean = Base.has as any;
    export const indexOf: (set: OrderedSet, x: number) => number = Base.indexOf as any;
    export const getAt: (set: OrderedSet, i: number) => number = Base.getAt as any;

    export const min: (set: OrderedSet) => number = Base.minValue as any;
    export const max: (set: OrderedSet) => number = Base.maxValue as any;
    export const size: (set: OrderedSet) => number = Base.size as any;
    export const hashCode: (set: OrderedSet) => number = Base.hashCode as any;

    export const areEqual: (a: OrderedSet, b: OrderedSet) => boolean = Base.areEqual as any;
    export const areIntersecting: (a: OrderedSet, b: OrderedSet) => boolean = Base.areIntersecting as any;
    export const isSubset: (a: OrderedSet, b: OrderedSet) => boolean = Base.isSubset as any;

    export const union: (a: OrderedSet, b: OrderedSet) => OrderedSet = Base.union as any;
    export const intersect: (a: OrderedSet, b: OrderedSet) => OrderedSet = Base.intersect as any;
    export const subtract: (a: OrderedSet, b: OrderedSet) => OrderedSet = Base.subtract as any;

    export const getPredIndex: (set: OrderedSet, x: number) => number = Base.getPredIndex as any;
    export const getPredIndex1: (set: OrderedSet, x: number, start: number, end: number) => number = Base.getPredIndex1 as any;
    export const getIntervalRange: (set: OrderedSet, min: number, max: number) => IndexRange = (set, min, max) => Base.getIntervalRange(set as any, min, max, IndexRange());
    export const getIntervalRange1: (set: OrderedSet, min: number, max: number, target: IndexRange) => IndexRange = Base.getIntervalRange as any;

    export const segments = SegmentIterator
}

interface OrderedSet { '@type': 'int-ordered-set' }

export default OrderedSet