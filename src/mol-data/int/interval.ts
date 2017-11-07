/**
 * Copyright (c) 2017 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as Impl from './impl/interval'

namespace Interval {
    export const Empty: Interval = Impl.Empty as any;

    export const ofSingleton: (value: number) => Interval = (v) => Impl.ofRange(v, v) as any;
    /** Create interval [min, max] */
    export const ofRange: (min: number, max: number) => Interval = Impl.ofRange as any;
    /** Create interval [min, max) */
    export const ofBounds: (start: number, end: number) => Interval = Impl.ofBounds as any;
    export const is: (v: any) => v is Interval = Impl.is as any;

    export const has: (interval: Interval, x: number) => boolean = Impl.has as any;
    export const indexOf: (interval: Interval, x: number) => number = Impl.indexOf as any;
    export const getAt: (interval: Interval, i: number) => number = Impl.getAt as any;

    export const start: (interval: Interval) => number = Impl.start as any;
    export const end: (interval: Interval) => number = Impl.end as any;
    export const min: (interval: Interval) => number = Impl.min as any;
    export const max: (interval: Interval) => number = Impl.max as any;
    export const size: (interval: Interval) => number = Impl.size as any;
    export const hashCode: (interval: Interval) => number = Impl.hashCode as any;

    export const areEqual: (a: Interval, b: Interval) => boolean = Impl.areEqual as any;
    export const areIntersecting: (a: Interval, b: Interval) => boolean = Impl.areIntersecting as any;
    export const isSubInterval: (a: Interval, b: Interval) => boolean = Impl.isSubInterval as any;

    export const findPredecessorIndex: (interval: Interval, x: number) => number = Impl.findPredecessorIndex as any;
    export const findPredecessorIndexInInterval: (interval: Interval, x: number, bounds: Interval) => number = Impl.findPredecessorIndexInInterval as any;
    export const findRange: (interval: Interval, min: number, max: number) => Interval = Impl.findRange as any;

    export const intersect: (a: Interval, b: Interval) => Interval = Impl.intersect as any;
}

interface Interval { '@type': 'int-interval' }

export default Interval