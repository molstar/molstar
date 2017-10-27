/**
 * Copyright (c) 2017 molio contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import IntTuple from './int-tuple'

/** Closed/Open inteval [a, b) (iterate as a <= i < b) */
namespace Interval {
    export const Empty: Interval = Impl.Empty as any;

    // export const create: (min: number, max: number) => Interval = Impl.create as any;

    // export const indexOf: (set: Interval, x: number) => number = Impl.indexOf as any;
    // export const getAt: (set: Interval, i: number) => number = Impl.getAt as any;

    // export const start: (set: Interval) => number = Impl.start as any;
    // export const end: (set: Interval) => number = Impl.end as any;

    // export const min: (set: Interval) => number = Impl.min as any;
    // export const max: (set: Interval) => number = Impl.max as any;
    // export const size: (set: Interval) => number = Impl.size as any;
    // export const hashCode: (set: Interval) => number = Impl.hashCode as any;

    // export const areEqual: (a: Interval, b: Interval) => boolean = Impl.areEqual as any;
    // export const areIntersecting: (a: Interval, b: Interval) => boolean = Impl.areIntersecting as any;
    // export const isSubset: (a: Interval, b: Interval) => boolean = Impl.isSubset as any;

}

interface Interval { '@type': 'int-interval' }

export default Interval

namespace Impl {export const Empty = IntTuple.Zero;

    export function create(min: number, max: number) { return max < min ? Empty : IntTuple.create(min, max + 1); }

    export const start = IntTuple.fst
    export const end = IntTuple.snd
    export const min = IntTuple.fst
    export function max(i: IntTuple) { return IntTuple.snd(i) + 1; }
    export function size(i: IntTuple) { return IntTuple.snd(i) - IntTuple.fst(i); }

    export function has(int: IntTuple, v: number) { return IntTuple.fst(int) <= v && v < IntTuple.snd(int); }
}