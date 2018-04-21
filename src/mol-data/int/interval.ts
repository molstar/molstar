/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/interval'

namespace Interval {
    export const Empty: Interval = Impl.Empty as any;

    export const ofSingleton: (value: number) => Interval = (v) => Impl.ofRange(v, v) as any;
    /** Create interval from range [min, max] */
    export const ofRange: (min: number, max: number) => Interval = Impl.ofRange as any;
    /** Create interval from bounds [start, end), i.e. [start, end - 1] */
    export const ofBounds: (start: number, end: number) => Interval = Impl.ofBounds as any;
    export const is: (v: any) => v is Interval = Impl.is as any;

    /** Test if a value is within the bounds of the interval */
    export const has: (interval: Interval, x: number) => boolean = Impl.has as any;
    export const indexOf: (interval: Interval, x: number) => number = Impl.indexOf as any;
    export const getAt: (interval: Interval, i: number) => number = Impl.getAt as any;

    /** Start value of the interval, same as min value */
    export const start: (interval: Interval) => number = Impl.start as any;
    /** End value of the interval, same as max + 1 */
    export const end: (interval: Interval) => number = Impl.end as any;
    /** Min value of the interval, same as start value */
    export const min: (interval: Interval) => number = Impl.min as any;
    /** Max value of the interval, same as end - 1 */
    export const max: (interval: Interval) => number = Impl.max as any;
    /** Number of values in the interval */
    export const size: (interval: Interval) => number = Impl.size as any;
    /** Hash code describing the interval */
    export const hashCode: (interval: Interval) => number = Impl.hashCode as any;

    /** Test if two intervals are identical */
    export const areEqual: (a: Interval, b: Interval) => boolean = Impl.areEqual as any;
    /** Test if two intervals are intersecting, i.e. their bounds overlap */
    export const areIntersecting: (a: Interval, b: Interval) => boolean = Impl.areIntersecting as any;

    /** Test if interval b is fully included in interval a */
    export const isSubInterval: (a: Interval, b: Interval) => boolean = Impl.isSubInterval as any;

    export const findPredecessorIndex: (interval: Interval, x: number) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: (interval: Interval, x: number, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: (interval: Interval, min: number, max: number) => Interval = Impl.findRange as any;

    /** Get a new interval that is the intersection of the two intervals */
    export const intersect: (a: Interval, b: Interval) => Interval = Impl.intersect as any;
}

/** Interval describing a range [min, max] of values */
interface Interval { '@type': 'int-interval' }

export default Interval