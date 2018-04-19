/*
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Result/*, Lookup3D*/ } from './common'
import { Box3D } from '../primitives/box3d';
import { Sphere3D } from '../primitives/sphere3d';
import { PositionData } from '../common';
import { Vec3 } from '../../linear-algebra';

// class GridLookup3D implements Lookup3D {
//     private result = Result.create();
//     private pivot = [0.1, 0.1, 0.1];
//     private radiusSq = 0.1;

//     find(x: number, y: number, z: number, radius: number): Result {
//         throw new Error("Method not implemented.");
//     }
//     check(x: number, y: number, z: number, radius: number): boolean {
//         throw new Error("Method not implemented.");
//     }
//     boundingBox: Box3D;
//     boundingSphere: Sphere3D;
// }

/**
 * Adapted from https://github.com/arose/ngl
 * MIT License Copyright (C) 2014+ Alexander Rose
 */

interface InputData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    radius?: ArrayLike<number>,
    indices: ArrayLike<number>
}

interface SpatialHash {
    size: number[],
    min: number[],
    grid: Uint32Array,
    delta: number[],
    bucketOffset: Uint32Array,
    bucketCounts: Int32Array,
    bucketArray: Int32Array,
    data: InputData,
    maxRadius: number,
    expandedBox: Box3D,
    boundingBox: Box3D,
    boundingSphere: Sphere3D
}

interface BuildState {
    size: number[],
    delta: number[],
    data: InputData,
    expandedBox: Box3D,
    boundingBox: Box3D,
    boundingSphere: Sphere3D,
    count: number
}

function _build(state: BuildState): SpatialHash {
    const { expandedBox, size: [sX, sY, sZ], data: { x: px, y: py, z: pz, radius, indices }, count, delta } = state;
    const n = sX * sY * sZ;
    const { min: [minX, minY, minZ] } = expandedBox;

    let maxRadius = 0;
    let bucketCount = 0;
    const grid = new Uint32Array(n);
    const bucketIndex = new Int32Array(count);
    for (let t = 0, _t = indices.length; t < _t; t++) {
        const i = indices[t];
        const x = Math.floor((px[i] - minX) / delta[0]);
        const y = Math.floor((py[i] - minY) / delta[1]);
        const z = Math.floor((pz[i] - minZ) / delta[2]);
        const idx = (((x * sY) + y) * sZ) + z;

        if ((grid[idx] += 1) === 1) {
            bucketCount += 1
        }
        bucketIndex[t] = idx;
    }

    if (radius) {
        for (let t = 0, _t = indices.length; t < _t; t++) {
            const i = indices[t];
            if (radius[i] > maxRadius) maxRadius = radius[i];
        }
    }

    const bucketCounts = new Int32Array(bucketCount);
    for (let i = 0, j = 0; i < n; i++) {
        const c = grid[i];
        if (c > 0) {
            grid[i] = j + 1;
            bucketCounts[j] = c;
            j += 1;
        }
    }

    const bucketOffset = new Uint32Array(count);
    for (let i = 1; i < count; ++i) {
        bucketOffset[i] += bucketOffset[i - 1] + bucketCounts[i - 1];
    }

    const bucketFill = new Int32Array(bucketCount);
    const bucketArray = new Int32Array(count);
    for (let i = 0; i < count; i++) {
        const bucketIdx = grid[bucketIndex[i]]
        if (bucketIdx > 0) {
            const k = bucketIdx - 1;
            bucketArray[bucketOffset[k] + bucketFill[k]] = i;
            bucketFill[k] += 1;
        }
    }

    return {
        size: state.size,
        bucketArray,
        bucketCounts,
        bucketOffset,
        grid,
        delta,
        min: state.expandedBox.min,
        data: state.data,
        maxRadius,
        expandedBox: state.expandedBox,
        boundingBox: state.boundingBox,
        boundingSphere: state.boundingSphere
    }
}

export function build(data: PositionData) {
    const boundingBox = Box3D.computeBounding(data);
    const expandedBox = Box3D.expand(boundingBox, Vec3.create(0.5, 0.5, 0.5));
    const boundingSphere = Sphere3D.computeBounding(data);

    // TODO: specialized code that does not require the indices annotation?
    const indices = data.indices ? data.indices as number[] : new Int32Array(data.x.length) as any as number[];
    if (!data.indices) {
        for (let i = 0, _i = data.x.length; i < _i; i++) indices[i] = i;
    }

    const S = Vec3.sub(Vec3.zero(), expandedBox.max, expandedBox.min);
    let delta, size;

    if (indices.length > 0) {
        // size of the box
        // required "grid volume" so that each cell contains on average 32 elements.
        const V = Math.ceil(indices.length / 32);
        const f = Math.pow(V / (S[0] * S[1] * S[2]), 1 / 3);
        size = [Math.ceil(S[0] * f), Math.ceil(S[1] * f), Math.ceil(S[2] * f)];
        delta = [S[0] / size[0], S[1] / size[1], S[2] / size[2]];
    } else {
        delta = S;
        size = [1, 1, 1];
    }

    const inputData: InputData = {
        x: data.x,
        y: data.y,
        z: data.z,
        indices,
        radius: data.radius
    };

    const state: BuildState = {
        size,
        data: inputData,
        expandedBox,
        boundingBox,
        boundingSphere,
        count: indices.length,
        delta
    }

    return _build(state);
}

export function inSphere(structure: SpatialHash, x: number, y: number, z: number, radius: number): Result {
    const { min, size: [sX, sY, sZ], bucketOffset, bucketCounts, bucketArray, grid, data: { x: px, y: py, z: pz, indices }, delta } = structure;
    //const { radius: r, radiusSq: rSq, pivot: [x, y, z], isCheck, mask } = ctx;
    const r = radius, rSq = r * r;

    const result = Result.create();

    const loX = Math.max(0, Math.floor((x - r - min[0]) / delta[0]));
    const loY = Math.max(0, Math.floor((y - r - min[1]) / delta[1]));
    const loZ = Math.max(0, Math.floor((z - r - min[2]) / delta[2]));

    const hiX = Math.min(sX - 1, Math.floor((x + r - min[0]) / delta[0]));
    const hiY = Math.min(sY - 1, Math.floor((y + r - min[1]) / delta[1]));
    const hiZ = Math.min(sZ - 1, Math.floor((z + r - min[2]) / delta[2]));

    if (loX > hiX || loY > hiY || loZ > hiZ) return result;

    for (let ix = loX; ix <= hiX; ix++) {
        for (let iy = loY; iy <= hiY; iy++) {
            for (let iz = loZ; iz <= hiZ; iz++) {
                const bucketIdx = grid[(((ix * sY) + iy) * sZ) + iz];
                if (bucketIdx === 0) continue;

                const k = bucketIdx - 1;
                const offset = bucketOffset[k];
                const count = bucketCounts[k];
                const end = offset + count;

                for (let i = offset; i < end; i++) {
                    const idx = indices[bucketArray[i]];

                    const dx = px[idx] - x;
                    const dy = py[idx] - y;
                    const dz = pz[idx] - z;
                    const distSq = dx * dx + dy * dy + dz * dz;

                    if (distSq <= rSq) {
                        // if (isCheck) return true;
                        Result.add(result, bucketArray[i], distSq);
                    }
                }
            }
        }
    }
    return result;
}