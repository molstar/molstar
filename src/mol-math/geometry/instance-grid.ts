/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OrderedSet } from '../../mol-data/int/ordered-set';
import { Box3D } from '../geometry';
import { Vec3 } from '../linear-algebra/3d/vec3';
import { PositionData } from './common';
import { GridLookup3D } from './lookup3d/grid';
import { Sphere3D } from './primitives/sphere3d';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4Offset = Vec3.transformMat4Offset;
const b3add = Box3D.add;

export type InstanceGrid = {
    readonly cellSize: number
    readonly cellCount: number
    readonly cellOffsets: Uint32Array
    readonly cellSpheres: Float32Array
    readonly cellTransform: Float32Array
    readonly cellInstance: Float32Array
}

export type InstanceData = {
    instanceCount: number
    instance: Float32Array
    transform: Float32Array
    invariantBoundingSphere: Sphere3D
}

export function createEmptyInstanceGrid(): InstanceGrid {
    return {
        cellSize: 0,
        cellCount: 0,
        cellOffsets: new Uint32Array(),
        cellSpheres: new Float32Array(),
        cellTransform: new Float32Array(),
        cellInstance: new Float32Array(),
    };
}

export function calcInstanceGrid(instanceData: InstanceData, cellSize: number): InstanceGrid {
    const { instanceCount, instance, transform, invariantBoundingSphere } = instanceData;
    // console.time('calcInstanceGrid grid');
    const x = new Float32Array(instanceCount);
    const y = new Float32Array(instanceCount);
    const z = new Float32Array(instanceCount);
    const indices = OrderedSet.ofBounds(0, instanceCount);

    const box = Box3D.setEmpty(Box3D());

    const { center, radius } = invariantBoundingSphere;
    const rv = Vec3.create(radius, radius, radius);

    const v = Vec3();
    for (let i = 0; i < instanceCount; ++i) {
        v3transformMat4Offset(v, center, transform, 0, 0, i * 16);
        x[i] = v[0];
        y[i] = v[1];
        z[i] = v[2];
        b3add(box, v);
    }
    Box3D.expand(box, box, rv);

    const positionData: PositionData = { x, y, z, indices };
    const boundary = { box, sphere: Sphere3D.fromBox3D(Sphere3D(), box) };
    const lookup = GridLookup3D(positionData, boundary, Vec3.create(cellSize, cellSize, cellSize));
    // console.timeEnd('calcInstanceGrid grid');

    const { array, offset, count } = lookup.buckets;

    const cellCount = offset.length;
    const cellOffsets = new Uint32Array(cellCount + 1);
    const cellSpheres = new Float32Array(cellCount * 4);
    const cellTransform = new Float32Array(instanceCount * 16);
    const cellInstance = new Float32Array(instanceCount);

    const b = Box3D();
    const s = Sphere3D();

    let k = 0;
    for (let i = 0; i < cellCount; ++i) {
        const start = offset[i];
        const size = count[i];
        cellOffsets[i] = start;
        const kStart = k;
        for (let j = start, jl = start + size; j < jl; ++j) {
            const idx = array[j];
            cellInstance[k] = instance[idx];
            for (let l = 0; l < 16; ++l) {
                cellTransform[k * 16 + l] = transform[idx * 16 + l];
            }
            k += 1;
        }

        if (size === 1) {
            v3transformMat4Offset(cellSpheres, center, cellTransform, i * 4, 0, kStart * 16);
            cellSpheres[i * 4 + 3] = radius;
        } else {
            Box3D.setEmpty(b);
            const o = kStart * 16;
            for (let l = 0; l < size; ++l) {
                v3transformMat4Offset(v, center, cellTransform, 0, 0, l * 16 + o);
                b3add(b, v);
            }
            Box3D.expand(b, b, rv);
            Sphere3D.fromBox3D(s, b);
            Sphere3D.toArray(s, cellSpheres, i * 4);
        }
    }
    cellOffsets[cellCount] = offset[cellCount - 1] + count[cellCount - 1];

    // console.log({
    //     cellSize,
    //     cellCount,
    //     cellOffsets,
    //     cellSpheres,
    //     cellTransform,
    //     cellInstance,
    // });

    return {
        cellSize,
        cellCount,
        cellOffsets,
        cellSpheres,
        cellTransform,
        cellInstance,
    };
}
