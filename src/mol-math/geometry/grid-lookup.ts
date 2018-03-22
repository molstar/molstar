/*
 * Copyright (c) 2017 MolQL contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// TODO: 3d grid lookup (use single cell for small sets), make bounding sphere part of the structure
// TODO: assign radius to points?

/**
 * A "masked" 3D spatial lookup structure.
 */

import Mask from 'mol-util/mask'

export type FindFunc = (x: number, y: number, z: number, radius: number) => Result
export type CheckFunc = (x: number, y: number, z: number, radius: number) => boolean

export interface Result {
    readonly count: number,
    readonly indices: number[],
    readonly squaredDistances: number[]
}

export interface Positions { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }

interface GridLookup { find: (mask: Mask) => FindFunc, check: (mask: Mask) => CheckFunc }

function GridLookup(positions: Positions): GridLookup {
    const tree = build(createInputData(positions));
    return {
        find(mask) {
            const ctx = QueryContext.create(tree, mask, false);
            return function (x: number, y: number, z: number, radius: number) {
                QueryContext.update(ctx, x, y, z, radius);
                nearest(ctx);
                return ctx.buffer;
            }
        },
        check(mask) {
            const ctx = QueryContext.create(tree, mask, true);
            return function (x: number, y: number, z: number, radius: number) {
                QueryContext.update(ctx, x, y, z, radius);
                return nearest(ctx);
            }
        }
    }
}

interface InputData {
    bounds: Box3D,
    count: number,
    positions: Positions
}

interface Box3D {
    min: number[],
    max: number[]
}

namespace Box3D {
    export function createInfinite(): Box3D {
        return {
            min: [Number.MAX_VALUE, Number.MAX_VALUE, Number.MAX_VALUE],
            max: [-Number.MAX_VALUE, -Number.MAX_VALUE, -Number.MAX_VALUE]
        }
    }
}

/**
* Query context. Handles the actual querying.
*/
interface QueryContext<T> {
    structure: T,
    pivot: number[],
    radius: number,
    radiusSq: number,
    buffer: QueryContext.Buffer,
    mask: Mask,
    isCheck: boolean
}

namespace QueryContext {
    export interface Buffer extends Result {
        count: number;
        indices: any[];
        squaredDistances: number[];
    }

    export function add<T>(ctx: QueryContext<T>, distSq: number, index: number) {
        const buffer = ctx.buffer;
        buffer.squaredDistances[buffer.count] = distSq;
        buffer.indices[buffer.count++] = index;
    }

    function resetBuffer(buffer: Buffer) { buffer.count = 0; }

    function createBuffer(): Buffer {
        return {
            indices: [],
            count: 0,
            squaredDistances: []
        }
    }

    /**
     * Query the tree and store the result to this.buffer. Overwrites the old result.
     */
    export function update<T>(ctx: QueryContext<T>, x: number, y: number, z: number, radius: number) {
        ctx.pivot[0] = x;
        ctx.pivot[1] = y;
        ctx.pivot[2] = z;
        ctx.radius = radius;
        ctx.radiusSq = radius * radius;
        resetBuffer(ctx.buffer);
    }

    export function create<T>(structure: T, mask: Mask, isCheck: boolean): QueryContext<T> {
        return {
            structure,
            buffer: createBuffer(),
            pivot: [0.1, 0.1, 0.1],
            radius: 1.1,
            radiusSq: 1.1 * 1.1,
            mask,
            isCheck
        }
    }
}

function createInputData(positions: { x: ArrayLike<number>, y: ArrayLike<number>, z: ArrayLike<number> }): InputData {
    const { x, y, z } = positions;
    const bounds = Box3D.createInfinite();
    const count = x.length;
    const { min, max } = bounds;

    for (let i = 0; i < count; i++) {
        min[0] = Math.min(x[i], min[0]);
        min[1] = Math.min(y[i], min[1]);
        min[2] = Math.min(z[i], min[2]);
        max[0] = Math.max(x[i], max[0]);
        max[1] = Math.max(y[i], max[1]);
        max[2] = Math.max(z[i], max[2]);
    }

    return { positions, bounds, count };
}

/**
 * Adapted from https://github.com/arose/ngl
 * MIT License Copyright (C) 2014+ Alexander Rose
 */

interface SpatialHash {
    size: number[],
    min: number[],
    grid: Uint32Array,
    bucketOffset: Uint32Array,
    bucketCounts: Int32Array,
    bucketArray: Int32Array,
    positions: Positions
}

interface State {
    size: number[],
    positions: Positions,
    bounds: Box3D,
    count: number
}

const enum Constants { Exp = 3 }

function nearest(ctx: QueryContext<SpatialHash>): boolean {
    const { min: [minX, minY, minZ], size: [sX, sY, sZ], bucketOffset, bucketCounts, bucketArray, grid, positions: { x: px, y: py, z: pz } } = ctx.structure;
    const { radius: r, radiusSq: rSq, pivot: [x, y, z], isCheck, mask } = ctx;

    const loX = Math.max(0, (x - r - minX) >> Constants.Exp);
    const loY = Math.max(0, (y - r - minY) >> Constants.Exp);
    const loZ = Math.max(0, (z - r - minZ) >> Constants.Exp);

    const hiX = Math.min(sX, (x + r - minX) >> Constants.Exp);
    const hiY = Math.min(sY, (y + r - minY) >> Constants.Exp);
    const hiZ = Math.min(sZ, (z + r - minZ) >> Constants.Exp);

    for (let ix = loX; ix <= hiX; ix++) {
        for (let iy = loY; iy <= hiY; iy++) {
            for (let iz = loZ; iz <= hiZ; iz++) {
                const bucketIdx = grid[(((ix * sY) + iy) * sZ) + iz];

                if (bucketIdx > 0) {
                    const k = bucketIdx - 1;
                    const offset = bucketOffset[k];
                    const count = bucketCounts[k];
                    const end = offset + count;

                    for (let i = offset; i < end; i++) {
                        const idx = bucketArray[i];
                        if (!mask.has(idx)) continue;

                        const dx = px[idx] - x;
                        const dy = py[idx] - y;
                        const dz = pz[idx] - z;
                        const distSq = dx * dx + dy * dy + dz * dz;

                        if (distSq <= rSq) {
                            if (isCheck) return true;
                            QueryContext.add(ctx, distSq, idx)
                        }
                    }
                }
            }
        }
    }
    return ctx.buffer.count > 0;
}

function _build(state: State): SpatialHash {
    const { bounds, size: [sX, sY, sZ], positions: { x: px, y: py, z: pz }, count } = state;
    const n = sX * sY * sZ;
    const { min: [minX, minY, minZ] } = bounds;

    let bucketCount = 0;
    const grid = new Uint32Array(n);
    const bucketIndex = new Int32Array(count);
    for (let i = 0; i < count; i++) {
        const x = (px[i] - minX) >> Constants.Exp;
        const y = (py[i] - minY) >> Constants.Exp;
        const z = (pz[i] - minZ) >> Constants.Exp;
        const idx = (((x * sY) + y) * sZ) + z;
        if ((grid[idx] += 1) === 1) {
            bucketCount += 1
        }
        bucketIndex[i] = idx;
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
        min: state.bounds.min,
        positions: state.positions
    }
}

function build({ positions, bounds, count }: InputData): SpatialHash {
    const size = [
        ((bounds.max[0] - bounds.min[0]) >> Constants.Exp) + 1,
        ((bounds.max[1] - bounds.min[1]) >> Constants.Exp) + 1,
        ((bounds.max[2] - bounds.min[2]) >> Constants.Exp) + 1
    ];

    const state: State = {
        size,
        positions,
        bounds,
        count
    }

    return _build(state);
}

export default GridLookup;