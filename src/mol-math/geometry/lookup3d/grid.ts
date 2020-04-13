/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Result, Lookup3D } from './common';
import { Box3D } from '../primitives/box3d';
import { Sphere3D } from '../primitives/sphere3d';
import { PositionData } from '../common';
import { Vec3 } from '../../linear-algebra';
import { OrderedSet } from '../../../mol-data/int';
import { Boundary } from '../../../mol-model/structure/structure/util/boundary';

interface GridLookup3D<T = number> extends Lookup3D<T> {
    readonly buckets: { readonly offset: ArrayLike<number>, readonly count: ArrayLike<number>, readonly array: ArrayLike<number> }
}

function GridLookup3D<T extends number = number>(data: PositionData, boundary: Boundary, cellSizeOrCount?: Vec3 | number): GridLookup3D<T> {
    return new GridLookup3DImpl<T>(data, boundary, cellSizeOrCount);
}

export { GridLookup3D };

class GridLookup3DImpl<T extends number = number> implements GridLookup3D<T> {
    private ctx: QueryContext<T>;
    boundary: Lookup3D['boundary'];
    buckets: GridLookup3D['buckets'];
    result: Result<T>

    find(x: number, y: number, z: number, radius: number): Result<T> {
        this.ctx.x = x;
        this.ctx.y = y;
        this.ctx.z = z;
        this.ctx.radius = radius;
        this.ctx.isCheck = false;
        query(this.ctx);
        return this.ctx.result;
    }

    check(x: number, y: number, z: number, radius: number): boolean {
        this.ctx.x = x;
        this.ctx.y = y;
        this.ctx.z = z;
        this.ctx.radius = radius;
        this.ctx.isCheck = true;
        return query(this.ctx);
    }

    constructor(data: PositionData, boundary: Boundary, cellSizeOrCount?: Vec3 | number) {
        const structure = build(data, boundary, cellSizeOrCount);
        this.ctx = createContext<T>(structure);
        this.boundary = { box: structure.boundingBox, sphere: structure.boundingSphere };
        this.buckets = { offset: structure.bucketOffset, count: structure.bucketCounts, array: structure.bucketArray };
        this.result = this.ctx.result;
    }
}

/**
 * Adapted from https://github.com/arose/ngl
 * MIT License Copyright (C) 2014+ Alexander Rose
 */

interface InputData {
    x: ArrayLike<number>,
    y: ArrayLike<number>,
    z: ArrayLike<number>,
    indices: OrderedSet,
    radius?: ArrayLike<number>,
}

interface Grid3D {
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
    elementCount: number
}

function _build(state: BuildState): Grid3D {
    const { expandedBox, size: [sX, sY, sZ], data: { x: px, y: py, z: pz, radius, indices }, elementCount, delta } = state;
    const n = sX * sY * sZ;
    const { min: [minX, minY, minZ] } = expandedBox;

    let maxRadius = 0;
    let bucketCount = 0;
    const grid = new Uint32Array(n);
    const bucketIndex = new Int32Array(elementCount);
    for (let t = 0; t < elementCount; t++) {
        const i = OrderedSet.getAt(indices, t);
        const x = Math.floor((px[i] - minX) / delta[0]);
        const y = Math.floor((py[i] - minY) / delta[1]);
        const z = Math.floor((pz[i] - minZ) / delta[2]);
        const idx = (((x * sY) + y) * sZ) + z;

        if ((grid[idx] += 1) === 1) {
            bucketCount += 1;
        }
        bucketIndex[t] = idx;
    }

    if (radius) {
        for (let t = 0; t < elementCount; t++) {
            const i = OrderedSet.getAt(indices, t);
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

    const bucketOffset = new Uint32Array(bucketCount);
    for (let i = 1; i < bucketCount; ++i) {
        bucketOffset[i] += bucketOffset[i - 1] + bucketCounts[i - 1];
    }

    const bucketFill = new Int32Array(bucketCount);
    const bucketArray = new Int32Array(elementCount);
    for (let i = 0; i < elementCount; i++) {
        const bucketIdx = grid[bucketIndex[i]];
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
    };
}

function build(data: PositionData, boundary: Boundary, cellSizeOrCount?: Vec3 | number) {
    // need to expand the grid bounds to avoid rounding errors
    const expandedBox = Box3D.expand(Box3D.empty(), boundary.box, Vec3.create(0.5, 0.5, 0.5));
    const { indices } = data;

    const S = Box3D.size(Vec3.zero(), expandedBox);
    let delta, size;

    const elementCount = OrderedSet.size(indices);

    const cellCount = typeof cellSizeOrCount === 'number' ? cellSizeOrCount : 32;
    const cellSize = Array.isArray(cellSizeOrCount) && cellSizeOrCount;

    if (cellSize) {
        size = [Math.ceil(S[0] / cellSize[0]), Math.ceil(S[1] / cellSize[1]), Math.ceil(S[2] / cellSize[2])];
        delta = cellSize;
    } else if (elementCount > 0) {
        // size of the box
        // required "grid volume" so that each cell contains on average 'cellCount' elements.
        const V = Math.ceil(elementCount / cellCount);
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
        boundingBox: boundary.box,
        boundingSphere: boundary.sphere,
        elementCount,
        delta
    };

    return _build(state);
}

interface QueryContext<T extends number = number> {
    grid: Grid3D,
    x: number,
    y: number,
    z: number,
    radius: number,
    result: Result<T>,
    isCheck: boolean
}

function createContext<T extends number = number>(grid: Grid3D): QueryContext<T> {
    return { grid, x: 0.1, y: 0.1, z: 0.1, radius: 0.1, result: Result.create(), isCheck: false };
}

function query<T extends number = number>(ctx: QueryContext<T>): boolean {
    const { min, size: [sX, sY, sZ], bucketOffset, bucketCounts, bucketArray, grid, data: { x: px, y: py, z: pz, indices, radius }, delta, maxRadius } = ctx.grid;
    const { radius: inputRadius, isCheck, x, y, z, result } = ctx;

    const r = inputRadius + maxRadius;
    const rSq = r * r;

    Result.reset(result);

    const loX = Math.max(0, Math.floor((x - r - min[0]) / delta[0]));
    const loY = Math.max(0, Math.floor((y - r - min[1]) / delta[1]));
    const loZ = Math.max(0, Math.floor((z - r - min[2]) / delta[2]));

    const hiX = Math.min(sX - 1, Math.floor((x + r - min[0]) / delta[0]));
    const hiY = Math.min(sY - 1, Math.floor((y + r - min[1]) / delta[1]));
    const hiZ = Math.min(sZ - 1, Math.floor((z + r - min[2]) / delta[2]));

    if (loX > hiX || loY > hiY || loZ > hiZ) return false;

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
                    const idx = OrderedSet.getAt(indices, bucketArray[i]);

                    const dx = px[idx] - x;
                    const dy = py[idx] - y;
                    const dz = pz[idx] - z;
                    const distSq = dx * dx + dy * dy + dz * dz;

                    if (distSq <= rSq) {
                        if (maxRadius > 0 && Math.sqrt(distSq) - radius![idx] > inputRadius) continue;
                        if (isCheck) return true;
                        Result.add(result, bucketArray[i], distSq);
                    }
                }
            }
        }
    }
    return result.count > 0;
}