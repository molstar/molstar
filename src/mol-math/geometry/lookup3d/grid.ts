/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { Result, Lookup3D } from './common';
import { Box3D } from '../primitives/box3d';
import { Sphere3D } from '../primitives/sphere3d';
import { PositionData } from '../common';
import { Vec3, EPSILON } from '../../linear-algebra';
import { OrderedSet } from '../../../mol-data/int';
import { Boundary } from '../boundary';
import { FibonacciHeap } from '../../../mol-util/fibonacci-heap';

interface GridLookup3D<T = number> extends Lookup3D<T> {
    readonly buckets: { readonly offset: ArrayLike<number>, readonly count: ArrayLike<number>, readonly array: ArrayLike<number> }
}

function GridLookup3D<T extends number = number>(data: PositionData, boundary: Boundary, cellSizeOrCount?: Vec3 | number): GridLookup3D<T> {
    return new GridLookup3DImpl<T>(data, boundary, cellSizeOrCount);
}

export { GridLookup3D };

class GridLookup3DImpl<T extends number = number> implements GridLookup3D<T> {
    private ctx: QueryContext;
    boundary: Lookup3D['boundary'];
    buckets: GridLookup3D['buckets'];
    result: Result<T>;

    find(x: number, y: number, z: number, radius: number, result?: Result<T>): Result<T> {
        this.ctx.x = x;
        this.ctx.y = y;
        this.ctx.z = z;
        this.ctx.radius = radius;
        this.ctx.isCheck = false;
        const ret = result ?? this.result;
        query(this.ctx, ret);
        return ret;
    }

    nearest(x: number, y: number, z: number, k: number = 1, result?: Result<T>): Result<T> {
        this.ctx.x = x;
        this.ctx.y = y;
        this.ctx.z = z;
        this.ctx.k = k;
        const ret = result ?? this.result;
        queryNearest(this.ctx, ret);
        return ret;
    }

    check(x: number, y: number, z: number, radius: number): boolean {
        this.ctx.x = x;
        this.ctx.y = y;
        this.ctx.z = z;
        this.ctx.radius = radius;
        this.ctx.isCheck = true;
        return query(this.ctx, this.result);
    }

    constructor(data: PositionData, boundary: Boundary, cellSizeOrCount?: Vec3 | number) {
        const structure = build(data, boundary, cellSizeOrCount);
        this.ctx = createContext(structure);
        this.boundary = { box: structure.boundingBox, sphere: structure.boundingSphere };
        this.buckets = { offset: structure.bucketOffset, count: structure.bucketCounts, array: structure.bucketArray };
        this.result = Result.create();
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
    const expandedBox = Box3D.expand(Box3D(), boundary.box, Vec3.create(0.5, 0.5, 0.5));
    const { indices } = data;

    const S = Box3D.size(Vec3(), expandedBox);
    let delta, size;

    const elementCount = OrderedSet.size(indices);

    const cellCount = typeof cellSizeOrCount === 'number' ? cellSizeOrCount : 32;
    const cellSize = Array.isArray(cellSizeOrCount) && cellSizeOrCount;

    if (cellSize && !Vec3.isZero(cellSize)) {
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

interface QueryContext {
    grid: Grid3D,
    x: number,
    y: number,
    z: number,
    k: number,
    radius: number,
    isCheck: boolean
}

function createContext(grid: Grid3D): QueryContext {
    return { grid, x: 0.1, y: 0.1, z: 0.1, k: 1, radius: 0.1, isCheck: false };
}

function query<T extends number = number>(ctx: QueryContext, result: Result<T>): boolean {
    const { min, size: [sX, sY, sZ], bucketOffset, bucketCounts, bucketArray, grid, data: { x: px, y: py, z: pz, indices, radius }, delta, maxRadius } = ctx.grid;
    const { radius: inputRadius, isCheck, x, y, z } = ctx;

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

const tmpDirVec = Vec3();
const tmpVec = Vec3();
const tmpMapG = new Map<number, boolean>();
const tmpArrG1 = [0.1];
const tmpArrG2 = [0.1];
const tmpHeapG = new FibonacciHeap();
function queryNearest<T extends number = number>(ctx: QueryContext, result: Result<T>): boolean {
    const { expandedBox: box, boundingSphere: { center }, size: [sX, sY, sZ], bucketOffset, bucketCounts, bucketArray, grid, data: { x: px, y: py, z: pz, indices, radius }, delta, maxRadius } = ctx.grid;
    const [minX, minY, minZ] = box.min;
    const { x, y, z, k } = ctx;
    const indicesCount = OrderedSet.size(indices);
    Result.reset(result);
    if (indicesCount === 0 || k <= 0) return false;
    let gX, gY, gZ;
    Vec3.set(tmpVec, x, y, z);
    if (!Box3D.containsVec3(box, tmpVec)) {
        // intersect ray pointing to box center
        Box3D.nearestIntersectionWithRay(tmpVec, box, tmpVec, Vec3.normalize(tmpDirVec, Vec3.sub(tmpDirVec, center, tmpVec)));
        gX = Math.max(0, Math.min(sX - 1, Math.floor((tmpVec[0] - minX) / delta[0])));
        gY = Math.max(0, Math.min(sY - 1, Math.floor((tmpVec[1] - minY) / delta[1])));
        gZ = Math.max(0, Math.min(sZ - 1, Math.floor((tmpVec[2] - minZ) / delta[2])));
    } else {
        gX = Math.floor((x - minX) / delta[0]);
        gY = Math.floor((y - minY) / delta[1]);
        gZ = Math.floor((z - minZ) / delta[2]);
    }
    let gCount = 1, nextGCount = 0, arrG = tmpArrG1, nextArrG = tmpArrG2, prevFurthestDist = Number.MAX_VALUE, prevNearestDist = -Number.MAX_VALUE, nearestDist = Number.MAX_VALUE, findFurthest = true, furthestDist = - Number.MAX_VALUE, distSqG: number;
    arrG.length = 0;
    nextArrG.length = 0;
    arrG.push(gX, gY, gZ);
    tmpMapG.clear();
    tmpHeapG.clear();
    while (result.count < indicesCount) {
        const arrGLen = gCount * 3;
        for (let ig = 0; ig < arrGLen; ig += 3) {
            gX = arrG[ig];
            gY = arrG[ig + 1];
            gZ = arrG[ig + 2];
            const gridId = (((gX * sY) + gY) * sZ) + gZ;
            if (tmpMapG.get(gridId)) continue; // already visited
            tmpMapG.set(gridId, true);
            distSqG = (gX - x) * (gX - x) + (gY - y) * (gY - y) + (gZ - z) * (gZ - z);
            if (!findFurthest && distSqG > furthestDist && distSqG < nearestDist) continue;
            // evaluate distances in the current grid point
            const bucketIdx = grid[gridId];
            if (bucketIdx !== 0) {
                const ki = bucketIdx - 1;
                const offset = bucketOffset[ki];
                const count = bucketCounts[ki];
                const end = offset + count;
                for (let i = offset; i < end; i++) {
                    const idx = OrderedSet.getAt(indices, bucketArray[i]);
                    const dx = px[idx] - x;
                    const dy = py[idx] - y;
                    const dz = pz[idx] - z;
                    let distSq = dx * dx + dy * dy + dz * dz;
                    if (maxRadius !== 0) {
                        const r = radius![idx];
                        const sqR = r * r;
                        if (findFurthest && distSq > furthestDist) furthestDist = distSq + sqR;
                        distSq = distSq - sqR;
                    } else {
                        if (findFurthest && distSq > furthestDist) furthestDist = distSq;
                    }
                    if (distSq > prevNearestDist && distSq <= furthestDist) {
                        tmpHeapG.insert(distSq, idx);
                        nearestDist = tmpHeapG.findMinimum()!.key as unknown as number;
                    }
                }
                if (prevFurthestDist < furthestDist) findFurthest = false;
            }
            // fill grid points array with valid adiacent positions
            for (let ix = -1; ix <= 1; ix++) {
                const xPos = gX + ix;
                if (xPos < 0 || xPos >= sX) continue;
                for (let iy = -1; iy <= 1; iy++) {
                    const yPos = gY + iy;
                    if (yPos < 0 || yPos >= sY) continue;
                    for (let iz = -1; iz <= 1; iz++) {
                        const zPos = gZ + iz;
                        if (zPos < 0 || zPos >= sZ) continue;
                        const gridId = (((xPos * sY) + yPos) * sZ) + zPos;
                        if (tmpMapG.get(gridId)) continue; // already visited
                        nextArrG.push(xPos, yPos, zPos);
                        nextGCount++;
                    }
                }
            }
        }
        if (nextGCount === 0) {
            while (!tmpHeapG.isEmpty() && result.count < k) {
                const node = tmpHeapG.extractMinimum();
                if (!node) throw new Error('Cannot extract minimum, should not happen');
                const { key: squaredDistance, value: index } = node;
                Result.add(result, index as number, squaredDistance as number);
            }
            if (result.count >= k) return result.count > 0;
            prevNearestDist = nearestDist;
            if (furthestDist === nearestDist) {
                findFurthest = true;
                prevFurthestDist = furthestDist;
                nearestDist = Number.MAX_VALUE;
            } else {
                nearestDist = furthestDist + EPSILON; // adding EPSILON fixes a bug
            }
            // resotre visibility of current gid points
            for (let ig = 0; ig < arrGLen; ig += 3) {
                gX = arrG[ig];
                gY = arrG[ig + 1];
                gZ = arrG[ig + 2];
                const gridId = (((gX * sY) + gY) * sZ) + gZ;
                tmpMapG.set(gridId, false);
            }
        } else {
            const tmp = arrG;
            arrG = nextArrG;
            nextArrG = tmp;
            nextArrG.length = 0;
            gCount = nextGCount;
            nextGCount = 0;
        }
    }
    return result.count > 0;
}
