/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as B from 'benchmark';
import { radixSort, sortArray } from '../mol-data/util/sort';

/**
 * Benchmarks the two grid-bucket compaction strategies used by `_build` in
 * src/mol-math/geometry/lookup3d/grid.ts, isolated from coordinate handling
 * (elements are given directly as precomputed cell ids):
 *
 * - scanAll: count elements per cell, then scan ALL n grid cells to assign
 *   bucket ids. Cost: O(elementCount + n).
 * - sortOccupied: additionally track occupied cell ids during the fill, then sort
 *   them and only touch occupied cells. Cost: O(elementCount + b log b) with
 *   b = number of occupied cells.
 *
 * Expectation: `sortOccupied` wins for sparse grids (b << n, e.g. large boxes hitting
 * the MaxVolume cap), `scanAll` wins for dense grids where b approaches n and
 * the O(b log b) sort exceeds the O(n) linear scan.
 *
 * `hybridStrategy` picks a sort path only for sparse grids and the scan path otherwise.
 * With radix sort the measured crossover is at ~25-30% occupancy => b <= n / 4.
 */

function genCellIds(n: number, elementCount: number, occupiedTarget: number) {
    const used = new Set<number>();
    while (used.size < occupiedTarget) used.add((Math.random() * n) | 0);
    const cells = new Int32Array(used.size);
    let i = 0;
    used.forEach(c => { cells[i++] = c; });

    const cellIds = new Int32Array(elementCount);
    // ensure every picked cell is occupied at least once
    for (let t = 0; t < cells.length; t++) cellIds[t] = cells[t];
    for (let t = cells.length; t < elementCount; t++) {
        cellIds[t] = cells[(Math.random() * cells.length) | 0];
    }
    return cellIds;
}

function scanAll(cellIds: Int32Array, n: number) {
    const elementCount = cellIds.length;
    let bucketCount = 0;
    const grid = new Uint32Array(n);
    for (let t = 0; t < elementCount; t++) {
        const idx = cellIds[t];
        if ((grid[idx] += 1) === 1) bucketCount += 1;
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
    return { grid, bucketCounts };
}

type Sorter = (occupied: Uint32Array, bucketCount: number, n: number) => ArrayLike<number>

const sortNative: Sorter = (occupied, bucketCount) => occupied.subarray(0, bucketCount).sort();

const sortQuick: Sorter = (occupied, bucketCount) => sortArray(occupied.subarray(0, bucketCount));

const sortRadix: Sorter = (occupied, bucketCount, n) => radixSort(occupied.subarray(0, bucketCount), n - 1);

function sortOccupied(cellIds: Int32Array, n: number, sorter: Sorter = sortNative) {
    const elementCount = cellIds.length;
    let bucketCount = 0;
    const grid = new Uint32Array(n);
    const occupied = new Uint32Array(Math.min(n, elementCount));
    for (let t = 0; t < elementCount; t++) {
        const idx = cellIds[t];
        if ((grid[idx] += 1) === 1) {
            occupied[bucketCount] = idx;
            bucketCount += 1;
        }
    }

    const occupiedSorted = sorter(occupied, bucketCount, n);
    const bucketCounts = new Int32Array(bucketCount);
    for (let j = 0; j < bucketCount; j++) {
        const i = occupiedSorted[j];
        bucketCounts[j] = grid[i];
        grid[i] = j + 1;
    }
    return { grid, bucketCounts };
}

function hybridStrategy(cellIds: Int32Array, n: number) {
    const elementCount = cellIds.length;
    let bucketCount = 0;
    const grid = new Uint32Array(n);
    const occupied = new Uint32Array(Math.min(n, elementCount));
    for (let t = 0; t < elementCount; t++) {
        const idx = cellIds[t];
        if ((grid[idx] += 1) === 1) {
            occupied[bucketCount] = idx;
            bucketCount += 1;
        }
    }

    const bucketCounts = new Int32Array(bucketCount);
    if (bucketCount <= n >>> 2) {
        const occupiedSorted = sortRadix(occupied, bucketCount, n);
        for (let j = 0; j < bucketCount; j++) {
            const i = occupiedSorted[j];
            bucketCounts[j] = grid[i];
            grid[i] = j + 1;
        }
    } else {
        for (let i = 0, j = 0; i < n; i++) {
            const c = grid[i];
            if (c > 0) {
                grid[i] = j + 1;
                bucketCounts[j] = c;
                j += 1;
            }
        }
    }
    return { grid, bucketCounts };
}

function verify(cellIds: Int32Array, n: number) {
    const a = scanAll(cellIds, n);
    const others = [
        sortOccupied(cellIds, n, sortNative),
        sortOccupied(cellIds, n, sortQuick),
        sortOccupied(cellIds, n, sortRadix),
        hybridStrategy(cellIds, n),
    ];
    for (const other of others) {
        if (a.bucketCounts.length !== other.bucketCounts.length) throw new Error('bucketCount mismatch');
        for (let i = 0; i < a.bucketCounts.length; i++) {
            if (a.bucketCounts[i] !== other.bucketCounts[i]) throw new Error(`bucketCounts mismatch at ${i}`);
        }
        for (let i = 0; i < n; i++) {
            if (a.grid[i] !== other.grid[i]) throw new Error(`grid mismatch at ${i}`);
        }
    }
}

function runTest(n: number, elementCount: number, occupiedTarget: number) {
    const cellIds = genCellIds(n, elementCount, occupiedTarget);
    verify(cellIds, n);

    const label = `n=${n} elements=${elementCount} occupied=${occupiedTarget} (${(100 * occupiedTarget / n).toFixed(2)}%)`;
    console.log(label);

    const suite = new B.Suite();
    suite
        .add('scan all cells', () => scanAll(cellIds, n))
        .add('sortOccupied (native .sort)', () => sortOccupied(cellIds, n, sortNative))
        .add('sortOccupied (sortArray quicksort)', () => sortOccupied(cellIds, n, sortQuick))
        .add('sortOccupied (radix sort)', () => sortOccupied(cellIds, n, sortRadix))
        .add('hybrid (radix if b <= n / 4)', () => hybridStrategy(cellIds, n))
        .on('cycle', (e: any) => console.log(`  ${String(e.target)}`))
        .on('complete', function (this: any) {
            const sorted = [this[0], this[1], this[2], this[3], this[4]].sort((a, b) => b.hz - a.hz);
            console.log(`  => fastest: ${sorted[0].name} (${(sorted[0].hz / sorted[1].hz).toFixed(2)}x over ${sorted[1].name})`);
        })
        .run();
    console.log('---------------------');
}

// sparse: large grid capped at MaxVolume, ~32 elements per occupied cell
// (typical big-structure case) => sortOccupied expected to win
runTest(2 ** 24, 100000, Math.floor(100000 / 32));
runTest(2 ** 24, 1000000, Math.floor(1000000 / 32));

// near the crossover (~1.5-2% occupancy)
runTest(2 ** 22, 200000, Math.floor(2 ** 22 * 0.015));

// dense: nearly every cell occupied, sort dominates => scanAll expected to win
runTest(2 ** 20, 1000000, Math.floor(2 ** 20 * 0.6));
runTest(2 ** 18, 500000, Math.floor(2 ** 18 * 0.9));

// tiny default grid (32^3): compaction negligible vs fill => roughly equal
runTest(32 ** 3, 100000, 32 ** 3 / 2);

// n=16777216 elements=100000 occupied=3125 (0.02%)
//   scan all cells x 28.40 ops/sec ±5.38% (49 runs sampled)
//   sortOccupied (native .sort) x 257 ops/sec ±1.53% (81 runs sampled)
//   sortOccupied (sortArray quicksort) x 270 ops/sec ±1.26% (86 runs sampled)
//   sortOccupied (radix sort) x 275 ops/sec ±1.24% (85 runs sampled)
//   hybrid (radix if b <= n / 4) x 261 ops/sec ±1.37% (82 runs sampled)
//   => fastest: sortOccupied (radix sort) (1.02x over sortOccupied (sortArray quicksort))
// ---------------------
// n=16777216 elements=1000000 occupied=31250 (0.19%)
//   scan all cells x 33.80 ops/sec ±0.95% (59 runs sampled)
//   sortOccupied (native .sort) x 52.17 ops/sec ±1.92% (65 runs sampled)
//   sortOccupied (sortArray quicksort) x 50.51 ops/sec ±1.76% (64 runs sampled)
//   sortOccupied (radix sort) x 53.30 ops/sec ±1.87% (55 runs sampled)
//   hybrid (radix if b <= n / 4) x 53.88 ops/sec ±2.35% (61 runs sampled)
//   => fastest: hybrid (radix if b <= n / 4) (1.01x over sortOccupied (radix sort))
// ---------------------
// n=4194304 elements=200000 occupied=62914 (1.50%)
//   scan all cells x 125 ops/sec ±0.76% (81 runs sampled)
//   sortOccupied (native .sort) x 122 ops/sec ±1.07% (80 runs sampled)
//   sortOccupied (sortArray quicksort) x 116 ops/sec ±0.72% (75 runs sampled)
//   sortOccupied (radix sort) x 171 ops/sec ±1.05% (80 runs sampled)
//   hybrid (radix if b <= n / 4) x 163 ops/sec ±1.26% (80 runs sampled)
//   => fastest: sortOccupied (radix sort) (1.04x over hybrid (radix if b <= n / 4))
// ---------------------
// n=1048576 elements=1000000 occupied=629145 (60.00%)
//   scan all cells x 127 ops/sec ±0.44% (82 runs sampled)
//   sortOccupied (native .sort) x 26.55 ops/sec ±0.47% (48 runs sampled)
//   sortOccupied (sortArray quicksort) x 23.39 ops/sec ±0.34% (43 runs sampled)
//   sortOccupied (radix sort) x 91.07 ops/sec ±0.78% (79 runs sampled)
//   hybrid (radix if b <= n / 4) x 109 ops/sec ±1.12% (81 runs sampled)
//   => fastest: scan all cells (1.17x over hybrid (radix if b <= n / 4))
// ---------------------
// n=262144 elements=500000 occupied=235929 (90.00%)
//   scan all cells x 516 ops/sec ±0.56% (92 runs sampled)
//   sortOccupied (native .sort) x 77.55 ops/sec ±0.35% (81 runs sampled)
//   sortOccupied (sortArray quicksort) x 67.41 ops/sec ±0.31% (71 runs sampled)
//   sortOccupied (radix sort) x 294 ops/sec ±0.45% (89 runs sampled)
//   hybrid (radix if b <= n / 4) x 472 ops/sec ±0.47% (90 runs sampled)
//   => fastest: scan all cells (1.09x over hybrid (radix if b <= n / 4))
// ---------------------
// n=32768 elements=100000 occupied=16384 (50.00%)
//   scan all cells x 4,837 ops/sec ±2.11% (82 runs sampled)
//   sortOccupied (native .sort) x 1,308 ops/sec ±1.05% (94 runs sampled)
//   sortOccupied (sortArray quicksort) x 1,150 ops/sec ±0.71% (97 runs sampled)
//   sortOccupied (radix sort) x 4,148 ops/sec ±2.18% (90 runs sampled)
//   hybrid (radix if b <= n / 4) x 5,230 ops/sec ±1.75% (91 runs sampled)
//   => fastest: hybrid (radix if b <= n / 4) (1.08x over scan all cells)
// ---------------------
