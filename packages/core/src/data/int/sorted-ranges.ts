/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Iterator as _Iterator } from '../iterator.js';
import { OrderedSet } from './ordered-set.js';
import type { Segmentation } from './segmentation.js';
import { SortedArray } from './sorted-array.js';

/** Pairs of min and max indices of sorted, non-overlapping ranges */
type SortedRanges<T extends number = number> = SortedArray<T>;

namespace SortedRanges {
  export function ofSortedRanges<T extends number = number>(array: ArrayLike<T>) {
    return SortedArray.ofSortedArray<T>(array);
  }
  export function start<T extends number = number>(ranges: SortedRanges<T>) {
    return ranges[0];
  }
  export function end<T extends number = number>(ranges: SortedRanges<T>) {
    return ranges[ranges.length - 1] + 1;
  }
  export function min<T extends number = number>(ranges: SortedRanges<T>) {
    return ranges[0];
  }
  export function max<T extends number = number>(ranges: SortedRanges<T>) {
    return ranges[ranges.length - 1];
  }
  export function size<T extends number = number>(ranges: SortedRanges<T>) {
    let size = 0;
    for (let i = 0, il = ranges.length; i < il; i += 2) {
      size += ranges[i + 1] - ranges[i] + 1;
    }
    return size;
  }
  export function count<T extends number = number>(ranges: SortedRanges<T>) {
    return ranges.length / 2;
  }

  export function startAt<T extends number = number>(ranges: SortedRanges<T>, index: number) {
    return ranges[index * 2];
  }
  export function endAt<T extends number = number>(ranges: SortedRanges<T>, index: number) {
    return ranges[index * 2 + 1] + 1;
  }

  export function minAt<T extends number = number>(ranges: SortedRanges<T>, index: number) {
    return ranges[index * 2];
  }
  export function maxAt<T extends number = number>(ranges: SortedRanges<T>, index: number) {
    return ranges[index * 2 + 1];
  }

  export function areEqual<T extends number = number>(a: SortedRanges<T>, b: SortedRanges<T>) {
    if (a.length !== b.length) return false;
    for (let i = 0, il = a.length; i < il; ++i) {
      if (a[i] !== b[i]) return false;
    }
    return true;
  }

  export function forEach<T extends number = number>(ranges: SortedRanges<T>, f: (value: T, i: number) => void) {
    let k = 0;
    for (let i = 0, il = ranges.length; i < il; i += 2) {
      for (let j = ranges[i], jl = ranges[i + 1]; j <= jl; ++j) {
        f(j, k);
        ++k;
      }
    }
  }

  /** Returns if a value of `set` is included in `ranges` */
  export function has<T extends number = number>(ranges: SortedRanges<T>, set: OrderedSet<T>) {
    return firstIntersectionIndex(ranges, set) !== -1;
  }

  /** Returns if a value of `set` is included in `ranges` from given index */
  export function hasFrom<T extends number = number>(ranges: SortedRanges<T>, set: OrderedSet<T>, from: number) {
    return firstIntersectionIndexFrom(ranges, set, from) !== -1;
  }

  export function firstIntersectionIndex<T extends number = number>(
    ranges: SortedRanges<T>,
    set: OrderedSet<T>,
  ): number {
    return firstIntersectionIndexFrom(ranges, set, 0);
  }

  export function firstIntersectionIndexFrom<T extends number = number>(
    ranges: SortedRanges<T>,
    set: OrderedSet<T>,
    from: number,
  ): number {
    const setMin = OrderedSet.min(set),
      setMax = OrderedSet.max(set);
    if (minAt(ranges, from) > setMax || max(ranges) < setMin) return -1;

    // the first range from `from` that ends at or after the set's min, by binary search: not a scan from `from`, an
    // Interval per range (perf-audit-2 01/P7)
    const il = count(ranges);
    let lo = from,
      hi = il;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (maxAt(ranges, mid) < setMin) lo = mid + 1;
      else hi = mid;
    }

    const setSize = OrderedSet.size(set);
    for (let i = lo; i < il; ++i) {
      const rMin = minAt(ranges, i);
      if (rMin > setMax) return -1;
      // the set's first value at or after the range's min, if within the range
      const idx = OrderedSet.findPredecessorIndex(set, rMin);
      if (idx < setSize && OrderedSet.getAt(set, idx) <= maxAt(ranges, i)) return i;
    }

    return -1;
  }

  export function transientSegments<T extends number = number, I extends number = number>(
    ranges: SortedRanges<T>,
    set: OrderedSet<T>,
  ) {
    return new Iterator<T, I>(ranges, set);
  }

  export class Iterator<T extends number = number, I extends number = number>
    implements _Iterator<Segmentation.Segment<I>>
  {
    private value: Segmentation.Segment<I> = { index: 0 as I, start: 0 as T, end: 0 as T };

    private curIndex = 0;

    hasNext: boolean = false;

    private updateValue() {
      this.value.index = this.curIndex as I;
      this.value.start = OrderedSet.findPredecessorIndex(this.set, startAt(this.ranges, this.curIndex));
      this.value.end = OrderedSet.findPredecessorIndex(this.set, endAt(this.ranges, this.curIndex));
    }

    move() {
      if (this.hasNext) {
        this.updateValue();
        this.curIndex = firstIntersectionIndexFrom(this.ranges, this.set, this.curIndex + 1);
        this.hasNext = this.curIndex !== -1;
      }
      return this.value;
    }

    constructor(
      private ranges: SortedRanges<T>,
      private set: OrderedSet<T>,
    ) {
      this.curIndex = firstIntersectionIndex(ranges, set);
      this.hasNext = this.curIndex !== -1;
    }
  }
}

export { SortedRanges };
