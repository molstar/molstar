/**
 * Copyright (c) 2023-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SortedArray } from '../../../mol-data/int';
import { ElementIndex } from '../../../mol-model/structure';
import { arrayExtend, range } from '../../../mol-util/array';


/** Represents a collection of disjoint elements ranges in a model (atoms, spheres, or gaussians).
 * The number of ranges is `ElementRanges.count(ranges)`,
 * the i-th range covers elements `[ranges.from[i], ranges.to[i])`. */
export interface ElementRanges {
    from: ElementIndex[],
    to: ElementIndex[],
}

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
        const concat = ElementRanges.empty();
        for (const r of ranges) {
            if (r) {
                arrayExtend(concat.from, r.from);
                arrayExtend(concat.to, r.to);
            }
        }
        const indices = range(concat.from.length).sort((i, j) => concat.from[i] - concat.from[j]); // sort by start of range
        const result = ElementRanges.empty();
        let last = -1;
        for (const i of indices) {
            const from = concat.from[i];
            const to = concat.to[i];
            if (last >= 0 && from <= result.to[last]) {
                if (to > result.to[last]) {
                    result.to[last] = to;
                }
            } else {
                result.from.push(from);
                result.to.push(to);
                last++;
            }
        }
        return result;
    },

    /** Return a sorted subset of `elements` which lie in any of `ranges` (i.e. set intersection of `elements` and `ranges`).
     * If `out` is provided, use it to store the result (clear any old contents).
     * If `outFirstElementIndex` is provided, fill `outFirstElementIndex.value` with the index of the first selected element (if any). */
    selectElementsInRanges(elements: SortedArray<ElementIndex>, ranges: ElementRanges, out?: ElementIndex[], outFirstElementIndex: { value?: number } = {}): ElementIndex[] {
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
