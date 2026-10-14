/**
 * Copyright (c) 2023-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SortedArray } from '@molstar/core/data/int';
import type { ElementIndex } from '../../../model/indexing.js';
import { MaxRetainedBufferCapacity } from '@molstar/core/util/int-buffer';

/** Represents a collection of disjoint elements ranges in a model (atoms, spheres, or gaussians).
 * The number of ranges is `ElementRanges.count(ranges)`,
 * the i-th range covers elements `[ranges.from[i], ranges.to[i])`. */
export interface ElementRanges {
  from: ElementIndex[];
  to: ElementIndex[];
}

/** Growable buffer of element ranges (not necessarily disjoint or sorted when filled via `append`), for reuse across calls.
 * Unlike `ElementRanges` it keeps an explicit `count` (see `IntBuffer` for why), so only the first `count` entries are valid. */
export class ElementRangesBuffer {
  from = new Int32Array(256);
  to = new Int32Array(256);
  count = 0;

  clear() {
    this.count = 0;
    return this;
  }

  private ensureCapacity(capacity: number) {
    if (capacity <= this.from.length) return;
    const newCapacity = Math.max(capacity, this.from.length * 2);
    const from = new Int32Array(newCapacity);
    const to = new Int32Array(newCapacity);
    from.set(this.from.subarray(0, this.count));
    to.set(this.to.subarray(0, this.count));
    this.from = from;
    this.to = to;
  }

  /** Same as `ElementRanges.add`: the added range must start after the end of the last range (touching ranges get merged). */
  add(from: number, to: number) {
    const n = this.count;
    if (n > 0) {
      const lastTo = this.to[n - 1];
      if (from < lastTo) throw new Error('Overlapping ranges not allowed');
      if (from === lastTo) {
        this.to[n - 1] = to;
        return this;
      }
    }
    this.ensureCapacity(n + 1);
    this.from[n] = from;
    this.to[n] = to;
    this.count = n + 1;
    return this;
  }

  /** Append ranges from `other` without any checks or merging (use `ElementRanges.normalize` on the result). */
  append(other: { from: ArrayLike<number>; to: ArrayLike<number>; count: number } | ElementRanges) {
    const n = 'count' in other ? other.count : other.from.length;
    this.ensureCapacity(this.count + n);
    for (let i = 0; i < n; i++) {
      this.from[this.count + i] = other.from[i];
      this.to[this.count + i] = other.to[i];
    }
    this.count += n;
    return this;
  }

  /** Copy the ranges into a new `ElementRanges`. */
  toRanges(): ElementRanges {
    return {
      from: Array.from(this.from.subarray(0, this.count)) as ElementIndex[],
      to: Array.from(this.to.subarray(0, this.count)) as ElementIndex[],
    };
  }

  /** Free the memory if the capacity grew above `limit` (so one huge call does not pin memory forever). */
  trim(limit: number) {
    if (this.from.length > limit) {
      this.from = new Int32Array(256);
      this.to = new Int32Array(256);
      this.count = 0;
    }
  }
}

/** Shared buffers (JS is single-threaded and none of the functions below is reentrant). Never returned to callers. */
const UnionBag = new ElementRangesBuffer();
let SweepStarts = new Int32Array(0);
let SweepEnds = new Int32Array(0);

export const ElementRanges = {
  /** Return the number of disjoined ranges in a `ElementRanges` object */
  count(ranges: ElementRanges): number {
    return ranges.from.length;
  },

  /** Create new `ElementRanges` without any elements */
  empty(): ElementRanges {
    return { from: [], to: [] };
  },

  /** Create new `ElementRanges` containing a single range of elements `[from, to)` */
  single(from: ElementIndex, to: ElementIndex): ElementRanges {
    return { from: [from], to: [to] };
  },

  /** Add a range of elements `[from, to)` to existing `ElementRanges` and return the modified original.
   * The added range must start after the end of the last existing range
   * (if it starts just on the next element, these two ranges will get merged). */
  add(ranges: ElementRanges, from: ElementIndex, to: ElementIndex): ElementRanges {
    const n = ElementRanges.count(ranges);
    if (n > 0) {
      const lastTo = ranges.to[n - 1];
      if (from < lastTo) throw new Error('Overlapping ranges not allowed');
      if (from === lastTo) {
        ranges.to[n - 1] = to;
      } else {
        ranges.from.push(from);
        ranges.to.push(to);
      }
    } else {
      ranges.from.push(from);
      ranges.to.push(to);
    }
    return ranges;
  },

  /** Apply function `func` to each range in `ranges` */
  foreach(ranges: ElementRanges, func: (from: ElementIndex, to: ElementIndex) => any) {
    const n = ElementRanges.count(ranges);
    for (let i = 0; i < n; i++) func(ranges.from[i], ranges.to[i]);
  },

  /** Apply function `func` to each range in `ranges` and return an array with results */
  map<T>(ranges: ElementRanges, func: (from: ElementIndex, to: ElementIndex) => T): T[] {
    const n = ElementRanges.count(ranges);
    const result: T[] = new Array(n);
    for (let i = 0; i < n; i++) result[i] = func(ranges.from[i], ranges.to[i]);
    return result;
  },

  /** Compute the set union of multiple `ElementRanges` objects (as sets of elements) */
  union(ranges: (ElementRanges | undefined)[]): ElementRanges {
    const bag = UnionBag.clear();
    for (const r of ranges) {
      if (r) bag.append(r);
    }
    const result = ElementRanges.normalize(bag);
    UnionBag.trim(MaxRetainedBufferCapacity);
    return result;
  },

  /** Turn a bag of ranges (in any order, possibly overlapping or touching) into a newly allocated `ElementRanges`
   * of disjoint ranges sorted by start (touching ranges are merged). `bag` is not modified.
   *
   * Instead of sorting ranges (which needs an index array and a comparator), the starts and the ends are sorted
   * independently as typed arrays (native, allocation-free numeric sort) and swept like interval-open/close events:
   * a merged range starts when the depth goes 0 -> 1 and ends when it goes 1 -> 0 (at equal coordinates a start is
   * processed before an end, so that touching ranges merge). */
  normalize(bag: ElementRangesBuffer): ElementRanges {
    const n = bag.count;
    const result = ElementRanges.empty();
    if (n === 0) return result;
    if (SweepStarts.length < n) {
      SweepStarts = new Int32Array(Math.max(1024, 1 << Math.ceil(Math.log2(n))));
      SweepEnds = new Int32Array(SweepStarts.length);
    }
    const starts = SweepStarts.subarray(0, n);
    const ends = SweepEnds.subarray(0, n);
    starts.set(bag.from.subarray(0, n));
    ends.set(bag.to.subarray(0, n));
    starts.sort();
    ends.sort();
    let i = 0;
    let j = 0;
    let depth = 0;
    while (j < n) {
      if (i < n && starts[i] <= ends[j]) {
        if (depth === 0) result.from.push(starts[i] as ElementIndex);
        depth++;
        i++;
      } else {
        depth--;
        if (depth === 0) result.to.push(ends[j] as ElementIndex);
        j++;
      }
    }
    if (SweepStarts.length > MaxRetainedBufferCapacity) {
      SweepStarts = new Int32Array(0);
      SweepEnds = new Int32Array(0);
    }
    return result;
  },

  /** Return a sorted subset of `elements` which lie in any of `ranges` (i.e. set intersection of `elements` and `ranges`).
   * If `out` is provided, use it to store the result (clear any old contents).
   * If `outFirstElementIndex` is provided, fill `outFirstElementIndex.value` with the index of the first selected element (if any). */
  selectElementsInRanges(
    elements: SortedArray<ElementIndex>,
    ranges: ElementRanges,
    out?: ElementIndex[],
    outFirstElementIndex: { value?: number } = {},
  ): ElementIndex[] {
    out ??= [];
    out.length = 0;
    outFirstElementIndex.value = undefined;

    const nElements = elements.length;
    const nRanges = ElementRanges.count(ranges);
    if (nElements <= nRanges) {
      // Implementation 1 (more efficient when there are fewer elements)
      let iRange = SortedArray.findPredecessorIndex(SortedArray.ofSortedArray(ranges.to), elements[0] + 1);
      for (let iElem = 0; iElem < nElements; iElem++) {
        const a = elements[iElem];
        while (iRange < nRanges && ranges.to[iRange] <= a) iRange++;
        const qualifies = iRange < nRanges && ranges.from[iRange] <= a;
        if (qualifies) {
          out.push(a);
          outFirstElementIndex.value ??= iElem;
        }
      }
    } else {
      // Implementation 2 (more efficient when there are fewer ranges)
      for (let iRange = 0; iRange < nRanges; iRange++) {
        const from = ranges.from[iRange];
        const to = ranges.to[iRange];
        for (let iElem = SortedArray.findPredecessorIndex(elements, from); iElem < nElements; iElem++) {
          const a = elements[iElem];
          if (a < to) {
            out.push(a);
            outFirstElementIndex.value ??= iElem;
          } else {
            break;
          }
        }
      }
    }
    return out;
  },
};
