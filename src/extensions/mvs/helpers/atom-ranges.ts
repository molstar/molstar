/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Adam Midlik <midlik@gmail.com>
 */

import { SortedArray } from '../../../mol-data/int';
import { ElementIndex } from '../../../mol-model/structure';
import { arrayExtend, range } from '../../../mol-util/array';


/** Represents a collection of disjoint atom ranges in a model.
 * The number of ranges is `AtomRanges.count(ranges)`,
 * the i-th range covers atoms `[ranges.from[i], ranges.to[i])`. */
export interface AtomRanges {
    from: ElementIndex[],
    to: ElementIndex[],
}

export const AtomRanges = {
    /** Return the number of disjoined ranges in a `AtomRanges` object */
    count(ranges: AtomRanges): number {
        return ranges.from.length;
    },

    /** Create new `AtomRanges` without any atoms */
    empty(): AtomRanges {
        return { from: [], to: [] };
    },

    /** Create new `AtomRanges` containing a single range of atoms `[from, to)` */
    single(from: ElementIndex, to: ElementIndex): AtomRanges {
        return { from: [from], to: [to] };
    },

    /** Add a range of atoms `[from, to)` to existing `AtomRanges` and return the modified original.
     * The added range must start after the end of the last existing range
     * (if it starts just on the next atom, these two ranges will get merged). */
    add(ranges: AtomRanges, from: ElementIndex, to: ElementIndex): AtomRanges {
        const n = AtomRanges.count(ranges);
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
    foreach(ranges: AtomRanges, func: (from: ElementIndex, to: ElementIndex) => any) {
        const n = AtomRanges.count(ranges);
        for (let i = 0; i < n; i++) func(ranges.from[i], ranges.to[i]);
    },

    /** Apply function `func` to each range in `ranges` and return an array with results */
    map<T>(ranges: AtomRanges, func: (from: ElementIndex, to: ElementIndex) => T): T[] {
        const n = AtomRanges.count(ranges);
        const result: T[] = new Array(n);
        for (let i = 0; i < n; i++) result[i] = func(ranges.from[i], ranges.to[i]);
        return result;
    },

    /** Compute the set union of multiple `AtomRanges` objects (as sets of atoms) */
    union(ranges: AtomRanges[]): AtomRanges {
        const concat = AtomRanges.empty();
        for (const r of ranges) {
            arrayExtend(concat.from, r.from);
            arrayExtend(concat.to, r.to);
        }
        const indices = range(concat.from.length).sort((i, j) => concat.from[i] - concat.from[j]); // sort by start of range
        const result = AtomRanges.empty();
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

    /** Return a sorted subset of `atoms` which lie in any of `ranges` (i.e. set intersection of `atoms` and `ranges`).
     * If `out` is provided, use it to store the result (clear any old contents).
     * If `outFirstAtomIndex` is provided, fill `outFirstAtomIndex.value` with the index of the first selected atom (if any). */
    selectAtomsInRanges(atoms: SortedArray<ElementIndex>, ranges: AtomRanges, out?: ElementIndex[], outFirstAtomIndex: { value?: number } = {}): ElementIndex[] {
        out ??= [];
        out.length = 0;
        outFirstAtomIndex.value = undefined;

        const nAtoms = atoms.length;
        const nRanges = AtomRanges.count(ranges);
        if (nAtoms <= nRanges) {
            // Implementation 1 (more efficient when there are fewer atoms)
            let iRange = SortedArray.findPredecessorIndex(SortedArray.ofSortedArray(ranges.to), atoms[0] + 1);
            for (let iAtom = 0; iAtom < nAtoms; iAtom++) {
                const a = atoms[iAtom];
                while (iRange < nRanges && ranges.to[iRange] <= a) iRange++;
                const qualifies = iRange < nRanges && ranges.from[iRange] <= a;
                if (qualifies) {
                    out.push(a);
                    outFirstAtomIndex.value ??= iAtom;
                }
            }
        } else {
            // Implementation 2 (more efficient when there are fewer ranges)
            for (let iRange = 0; iRange < nRanges; iRange++) {
                const from = ranges.from[iRange];
                const to = ranges.to[iRange];
                for (let iAtom = SortedArray.findPredecessorIndex(atoms, from); iAtom < nAtoms; iAtom++) {
                    const a = atoms[iAtom];
                    if (a < to) {
                        out.push(a);
                        outFirstAtomIndex.value ??= iAtom;
                    } else {
                        break;
                    }
                }
            }
        }
        return out;
    },
};
