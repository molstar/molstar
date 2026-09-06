/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../linear-algebra/3d/mat4';
import { Vec3 } from '../../linear-algebra/3d/vec3';
import { Sphere3D } from '../primitives/sphere3d';
import { calcInstanceGrid, createEmptyInstanceGrid, InstanceData, InstanceGrid } from '../instance-grid';
import { fillSerial } from '../../../mol-util/array';

function createInstanceData(positions: Vec3[], radius = 1): InstanceData {
    const instanceCount = positions.length;
    const transform = new Float32Array(instanceCount * 16);
    const m = Mat4();
    for (let i = 0; i < instanceCount; ++i) {
        Mat4.fromTranslation(m, positions[i]);
        transform.set(m, i * 16);
    }
    return {
        instanceCount,
        instance: fillSerial(new Float32Array(instanceCount)),
        transform,
        invariantBoundingSphere: Sphere3D.create(Vec3.create(0, 0, 0), radius),
    };
}

function checkInvariants(grid: InstanceGrid, data: InstanceData) {
    const { instanceCount, transform } = data;

    // offsets are monotonic and cover all instances
    expect(grid.cellOffsets.length).toBe(grid.cellCount + 1);
    expect(grid.cellOffsets[0]).toBe(0);
    for (let i = 0; i < grid.cellCount; ++i) {
        expect(grid.cellOffsets[i + 1]).toBeGreaterThan(grid.cellOffsets[i]);
    }
    expect(grid.cellOffsets[grid.cellCount]).toBe(instanceCount);

    // every instance appears exactly once
    const seen = new Set<number>();
    for (let i = 0; i < instanceCount; ++i) seen.add(grid.cellInstance[i]);
    expect(seen.size).toBe(instanceCount);

    // cellTransform is the original transform of each instance, in cell order
    for (let i = 0; i < instanceCount; ++i) {
        const id = grid.cellInstance[i];
        for (let k = 0; k < 16; ++k) {
            expect(grid.cellTransform[i * 16 + k]).toBe(transform[id * 16 + k]);
        }
    }

    // every instance center is contained in the sphere of its cell
    const { center, radius } = data.invariantBoundingSphere;
    const p = Vec3();
    for (let c = 0; c < grid.cellCount; ++c) {
        for (let i = grid.cellOffsets[c], il = grid.cellOffsets[c + 1]; i < il; ++i) {
            Vec3.transformMat4Offset(p, center, grid.cellTransform as unknown as number[], 0, 0, i * 16);
            const dx = p[0] - grid.cellSpheres[c * 4];
            const dy = p[1] - grid.cellSpheres[c * 4 + 1];
            const dz = p[2] - grid.cellSpheres[c * 4 + 2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            expect(dist).toBeLessThanOrEqual(grid.cellSpheres[c * 4 + 3] + 1e-4);
        }
    }

    // batches partition the cells; after reorder batchCell is serial
    expect(grid.batchOffsets.length).toBe(grid.batchCount + 1);
    expect(grid.batchOffsets[0]).toBe(0);
    expect(grid.batchOffsets[grid.batchCount]).toBe(grid.cellCount);
    expect(Array.from(grid.batchCell)).toEqual(Array.from(fillSerial(new Uint32Array(grid.cellCount))));

    // every cell sphere is contained in the sphere of its batch
    for (let b = 0; b < grid.batchCount; ++b) {
        for (let c = grid.batchOffsets[b], cl = grid.batchOffsets[b + 1]; c < cl; ++c) {
            const dx = grid.cellSpheres[c * 4] - grid.batchSpheres[b * 4];
            const dy = grid.cellSpheres[c * 4 + 1] - grid.batchSpheres[b * 4 + 1];
            const dz = grid.cellSpheres[c * 4 + 2] - grid.batchSpheres[b * 4 + 2];
            const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
            expect(dist).toBeLessThanOrEqual(grid.batchSpheres[b * 4 + 3] + 1e-4);
        }
    }

    // radius must be accounted for in cell spheres
    for (let c = 0; c < grid.cellCount; ++c) {
        expect(grid.cellSpheres[c * 4 + 3]).toBeGreaterThanOrEqual(radius - 1e-4);
    }
}

describe('instance-grid', () => {
    it('createEmptyInstanceGrid', () => {
        const grid = createEmptyInstanceGrid();
        expect(grid.cellCount).toBe(0);
        expect(grid.batchCount).toBe(0);
        expect(grid.cellOffsets.length).toBe(0);
    });

    it('single instance', () => {
        const data = createInstanceData([Vec3.create(1, 2, 3)]);
        const grid = calcInstanceGrid(data, 10, 100);

        expect(grid.cellCount).toBe(1);
        expect(grid.batchCount).toBe(1);
        checkInvariants(grid, data);

        // the cell sphere is the transformed invariant sphere
        expect(grid.cellSpheres[0]).toBeCloseTo(1);
        expect(grid.cellSpheres[1]).toBeCloseTo(2);
        expect(grid.cellSpheres[2]).toBeCloseTo(3);
        expect(grid.cellSpheres[3]).toBeCloseTo(1);
    });

    it('two clusters far apart get separate cells', () => {
        const data = createInstanceData([
            Vec3.create(0, 0, 0), Vec3.create(1, 0, 0),
            Vec3.create(100, 0, 0), Vec3.create(101, 0, 0),
        ]);
        const grid = calcInstanceGrid(data, 10, 1000);

        expect(grid.cellCount).toBeGreaterThanOrEqual(2);
        checkInvariants(grid, data);

        // instances of the same cluster share a cell
        const cellOf = new Map<number, number>();
        for (let c = 0; c < grid.cellCount; ++c) {
            for (let i = grid.cellOffsets[c], il = grid.cellOffsets[c + 1]; i < il; ++i) {
                cellOf.set(grid.cellInstance[i], c);
            }
        }
        expect(cellOf.get(0)).toBe(cellOf.get(1));
        expect(cellOf.get(2)).toBe(cellOf.get(3));
        expect(cellOf.get(0)).not.toBe(cellOf.get(2));
    });

    it('many random instances keep invariants', () => {
        // deterministic LCG so the test is reproducible
        let seed = 7;
        const rand = () => (seed = (seed * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff;

        const positions: Vec3[] = [];
        for (let i = 0; i < 500; ++i) {
            positions.push(Vec3.create(rand() * 200, rand() * 200, rand() * 200));
        }
        const data = createInstanceData(positions, 2);
        const grid = calcInstanceGrid(data, 25, 100);

        expect(grid.cellCount).toBeGreaterThan(1);
        expect(grid.batchCount).toBeGreaterThan(1);
        checkInvariants(grid, data);
    });
});
