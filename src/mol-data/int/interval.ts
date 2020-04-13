/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/interval';

namespace Interval {
    export const Empty: Interval = Impl.Empty as any;

    export const ofSingleton: <T extends number = number>(value: T) => Interval<T> = (v) => Impl.ofRange(v, v) as any;
    /** Create interval from range [min, max] */
    export const ofRange: <T extends number = number>(min: T, max: T) => Interval<T> = Impl.ofRange as any;
    /** Create interval from bounds [start, end), i.e. [start, end - 1] */
    export const ofBounds: <T extends number = number>(start: T, end: T) => Interval<T> = Impl.ofBounds as any;
    /** Create interval from length [0, length), i.e. [0, length - 1] */
    export const ofLength: <T extends number = number>(length: T) => Interval<T> = Impl.ofLength as any;
    export const is: <T extends number = number>(v: any) => v is Interval<T> = Impl.is as any;

    /** Test if a value is within the bounds of the interval */
    export const has: <T extends number = number>(interval: Interval<T>, x: T) => boolean = Impl.has as any;
    /** Returns the index of `x` in `set` or -1 if not found. */
    export const indexOf: <T extends number = number>(interval: Interval<T>, x: T) => number = Impl.indexOf as any;
    export const getAt: <T extends number = number>(interval: Interval<T>, i: number) => T = Impl.getAt as any;

    /** Start value of the Interval<T>, same as min value */
    export const start: <T extends number = number>(interval: Interval<T>) => T = Impl.start as any;
    /** End value of the Interval<T>, same as max + 1 */
    export const end: <T extends number = number>(interval: Interval<T>) => T = Impl.end as any;
    /** Min value of the Interval<T>, same as start value */
    export const min: <T extends number = number>(interval: Interval<T>) => T = Impl.min as any;
    /** Max value of the Interval<T>, same as end - 1 */
    export const max: <T extends number = number>(interval: Interval<T>) => T = Impl.max as any;
    /** Number of values in the interval */
    export const size: <T extends number = number>(interval: Interval<T>) => number = Impl.size as any;
    /** Hash code describing the interval */
    export const hashCode: <T extends number = number>(interval: Interval<T>) => number = Impl.hashCode as any;
    /** String representation of the interval */
    export const toString: <T extends number = number>(interval: Interval<T>) => string = Impl.toString as any;

    /** Test if two intervals are identical */
    export const areEqual: <T extends number = number>(a: Interval<T>, b: Interval<T>) => boolean = Impl.areEqual as any;
    /** Test if two intervals are intersecting, i.e. their bounds overlap */
    export const areIntersecting: <T extends number = number>(a: Interval<T>, b: Interval<T>) => boolean = Impl.areIntersecting as any;

    /** Test if interval b is fully included in interval a */
    export const isSubInterval: <T extends number = number>(a: Interval<T>, b: Interval<T>) => boolean = Impl.isSubInterval as any;

    export const findPredecessorIndex: <T extends number = number>(interval: Interval<T>, x: T) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: <T extends number = number>(interval: Interval<T>, x: T, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: <T extends number = number>(interval: Interval<T>, min: T, max: T) => Interval = Impl.findRange as any;
    /** Size of the intersection of the two intervals */
    export const intersectionSize: <T extends number = number>(a: Interval<T>, b: Interval<T>) => number = Impl.intersectionSize as any;

    /** Get a new interval that is the intersection of the two intervals */
    export const intersect: <T extends number = number>(a: Interval<T>, b: Interval<T>) => Interval<T> = Impl.intersect as any;
}

/** Interval describing a range [min, max] of values */
interface Interval<T extends number = number> { '@type': 'int-interval' }

export default Interval;