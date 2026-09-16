/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GridLookup3D } from '../../geometry';
import { sortArray } from '../../../mol-data/util/sort';
import { OrderedSet } from '../../../mol-data/int/ordered-set';
import { Vec3 } from '../../linear-algebra/3d/vec3';
import { getBoundary } from '../boundary';

const xs = [0, 0, 1];
const ys = [0, 1, 0];
const zs = [0, 0, 0];
const rs = [0, 0.5, 1 / 3];

describe('GridLookup3d', () => {
    it('basic', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        let r = grid.find(0, 0, 0, 0);
        expect(r.count).toBe(1);
        expect(r.indices[0]).toBe(0);

        r = grid.nearest(0, 0, 0, 1);
        expect(r.count).toBe(1);
        expect(r.indices[0]).toBe(0);

        r = grid.find(0, 0, 0, 1);
        expect(r.count).toBe(3);
        expect(sortArray(r.indices)).toEqual([0, 1, 2]);

        r = grid.nearest(0, 0, 0, 3);
        expect(r.count).toBe(3);
        expect(sortArray(r.indices)).toEqual([0, 1, 2]);
    });

    it('radius', () => {
        const position = { x: xs, y: ys, z: zs, radius: [0, 0.5, 1 / 3], indices: OrderedSet.ofBounds(0, 3) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        let r = grid.find(0, 0, 0, 0);
        expect(r.count).toBe(1);
        expect(r.indices[0]).toBe(0);

        r = grid.nearest(0, 0, 0, 1);
        expect(r.count).toBe(1);
        expect(r.indices[0]).toBe(0);

        r = grid.find(0, 0, 0, 0.5);
        expect(r.count).toBe(2);
        expect(sortArray(r.indices)).toEqual([0, 1]);

        r = grid.nearest(0, 0, 0, 3);
        expect(r.count).toBe(3);
        expect(sortArray(r.indices)).toEqual([0, 1, 2]);
    });

    it('indexed', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofSingleton(1), radius: rs };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        let r = grid.find(0, 0, 0, 0);
        expect(r.count).toBe(0);

        r = grid.nearest(0, 0, 0, 1);
        expect(r.count).toBe(1);

        r = grid.find(0, 0, 0, 0.5);
        expect(r.count).toBe(1);
        expect(sortArray(r.indices)).toEqual([0]);

        r = grid.nearest(0, 0, 0, 3);
        expect(r.count).toBe(1);
        expect(sortArray(r.indices)).toEqual([0]);
    });

    it('sorted-array indices', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofSortedArray([0, 2]), radius: rs };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        let r = grid.find(0, 0, 0, 0);
        expect(r.count).toBe(1);
        expect(r.indices[0]).toBe(0);

        // element 1 (y = 1) is excluded, so within radius 1 only elements 0 and 2 are found
        r = grid.find(0, 0, 0, 1);
        expect(r.count).toBe(2);
        expect(sortArray(r.indices)).toEqual([0, 1]);

        r = grid.nearest(0, 0, 0, 2);
        expect(r.count).toBe(2);
        expect(sortArray(r.indices)).toEqual([0, 1]);
    });

    it('check', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        expect(grid.check(0, 0, 0, 0.5)).toBe(true);
        expect(grid.check(5, 5, 5, 0.5)).toBe(false);
    });

    it('approxNearest', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        let r = grid.approxNearest(0, 0, 0, 0.5);
        expect(r.count).toBe(1);
        expect(r.indices[0]).toBe(0);

        r = grid.approxNearest(5, 5, 5, 0.5);
        expect(r.count).toBe(0);
    });

    it('empty', () => {
        const position = { x: [] as number[], y: [] as number[], z: [] as number[], indices: OrderedSet.ofBounds(0, 0) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        expect(grid.find(0, 0, 0, 1).count).toBe(0);
        expect(grid.nearest(0, 0, 0, 1).count).toBe(0);
        expect(grid.check(0, 0, 0, 1)).toBe(false);
    });

    it('nearest from outside the boundary', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        const r = grid.nearest(50, 50, 50, 1);
        expect(r.count).toBe(1);
        // element 1 (0,1,0) and element 2 (1,0,0) are equidistant candidates; (0,1,0) wins on distance ties by heap order
        expect(r.squaredDistances[0]).toBeCloseTo(Math.min(
            (50 - 0) ** 2 + (50 - 1) ** 2 + 50 ** 2,
            (50 - 1) ** 2 + 50 ** 2 + 50 ** 2,
        ));
    });

    it('nearest with stopIf', () => {
        const position = { x: xs, y: ys, z: zs, indices: OrderedSet.ofBounds(0, 3) };
        const boundary = getBoundary(position);
        const grid = GridLookup3D(position, boundary);

        const r = grid.nearest(0, 0, 0, 3, () => true);
        expect(r.count).toBeGreaterThanOrEqual(1);
        expect(r.count).toBeLessThanOrEqual(3);
    });

    it('many points with fixed cell size match brute force', () => {
        // deterministic LCG so the test is reproducible
        let seed = 42;
        const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;

        const count = 1000;
        const px = new Float32Array(count);
        const py = new Float32Array(count);
        const pz = new Float32Array(count);
        for (let i = 0; i < count; ++i) {
            px[i] = rand() * 100;
            py[i] = rand() * 100;
            pz[i] = rand() * 100;
        }

        const position = { x: px, y: py, z: pz, indices: OrderedSet.ofBounds(0, count) };
        const boundary = getBoundary(position);
        // sparse fixed-size cells exercise the occupied-cell compaction path
        const grid = GridLookup3D(position, boundary, Vec3.create(3, 3, 3));

        const queries: [number, number, number, number][] = [[50, 50, 50, 10], [0, 0, 0, 25], [100, 100, 100, 15], [50, 0, 100, 5]];
        for (const [qx, qy, qz, qr] of queries) {
            const expected: number[] = [];
            for (let i = 0; i < count; ++i) {
                const dx = px[i] - qx, dy = py[i] - qy, dz = pz[i] - qz;
                if (dx * dx + dy * dy + dz * dz <= qr * qr) expected.push(i);
            }
            const r = grid.find(qx, qy, qz, qr);
            expect(sortArray(r.indices.slice(0, r.count))).toEqual(expected);
            expect(grid.check(qx, qy, qz, qr)).toBe(expected.length > 0);
        }
    });
});
