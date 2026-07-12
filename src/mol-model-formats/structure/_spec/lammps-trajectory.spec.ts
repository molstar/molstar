/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Ludovic Autin <autin@scripps.edu>
 */

import { Column } from '../../../mol-data/db';
import { getCanonicalOrder } from '../lammps-trajectory';

function frame(atomIds: number[]) {
    return { count: atomIds.length, atomId: Column.ofIntArray(atomIds) };
}

describe('lammps-trajectory getCanonicalOrder', () => {
    it('returns undefined when rows are already in id order (identity fast-path)', () => {
        expect(getCanonicalOrder(frame([1, 2, 3, 4]))).toBeUndefined();
        expect(getCanonicalOrder(frame([1]))).toBeUndefined();
        expect(getCanonicalOrder(frame([]))).toBeUndefined(); // empty frame
    });

    it('maps each row to the 0-based slot id - 1', () => {
        // row 0 -> id 3 -> slot 2, row 1 -> id 1 -> slot 0, row 2 -> id 2 -> slot 1
        expect(getCanonicalOrder(frame([3, 1, 2]))).toEqual(new Int32Array([2, 0, 1]));
        // out of order from the first row
        expect(getCanonicalOrder(frame([2, 1]))).toEqual(new Int32Array([1, 0]));
    });

    it('preserves the leading identity prefix before the first out-of-order row', () => {
        // rows 0,1 already in order, row 2/3 swapped
        expect(getCanonicalOrder(frame([1, 2, 4, 3]))).toEqual(new Int32Array([0, 1, 3, 2]));
    });

    it('returns undefined for a duplicate id landing in the already-filled prefix (dst < j)', () => {
        // identity prefix [1,2], then id 1 points back into it
        expect(getCanonicalOrder(frame([1, 2, 1]))).toBeUndefined();
    });

    it('returns undefined for a duplicate id after the permutation is allocated (seen[dst])', () => {
        // first row already out of order (allocates), then a repeated id
        expect(getCanonicalOrder(frame([2, 1, 1]))).toBeUndefined();
    });

    it('returns undefined for out-of-range ids', () => {
        expect(getCanonicalOrder(frame([0, 1, 2]))).toBeUndefined(); // id 0 -> slot -1
        expect(getCanonicalOrder(frame([1, 2, 4]))).toBeUndefined(); // id 4 -> slot 3 >= count
    });
});
