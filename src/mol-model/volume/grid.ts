/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SpacegroupCell, Box3D, Sphere3D } from '../../mol-math/geometry';
import { Tensor, Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Histogram, calculateHistogram } from '../../mol-math/histogram';
import { lerp } from '../../mol-math/interpolate';

// avoiding namespace lookup improved performance in Chrome (Aug 2020)
const v3transformMat4 = Vec3.transformMat4;
const v3lerp = Vec3.lerp;

/** The basic unit cell that contains the grid data. */
interface Grid {
    readonly transform: Grid.Transform,
    readonly cells: Tensor,
    readonly stats: Readonly<{
        min: number,
        max: number,
        mean: number,
        sigma: number
    }>
}

namespace Grid {
    export const One: Grid = {
        transform: { kind: 'matrix', matrix: Mat4.identity() },
        cells: Tensor.create(Tensor.Space([1, 1, 1], [0, 1, 2]), Tensor.Data1([0])),
        stats: { min: 0, max: 0, mean: 0, sigma: 0 },
    };

    export type Transform = { kind: 'spacegroup', cell: SpacegroupCell, fractionalBox: Box3D } | { kind: 'matrix', matrix: Mat4 }

    const _scale = Mat4(), _translate = Mat4();
    export function getGridToCartesianTransform(grid: Grid) {
        if (grid.transform.kind === 'matrix') {
            return Mat4.copy(Mat4(), grid.transform.matrix);
        }

        if (grid.transform.kind === 'spacegroup') {
            const { cells: { space } } = grid;
            const scale = Mat4.fromScaling(_scale, Vec3.div(Vec3(), Box3D.size(Vec3(), grid.transform.fractionalBox), Vec3.ofArray(space.dimensions)));
            const translate = Mat4.fromTranslation(_translate, grid.transform.fractionalBox.min);
            return Mat4.mul3(Mat4(), grid.transform.cell.fromFractional, translate, scale);
        }

        return Mat4.identity();
    }

    export function areEquivalent(gridA: Grid, gridB: Grid) {
        return gridA === gridB;
    }

    export function isEmpty(grid: Grid) {
        return grid.cells.data.length === 0;
    }

    export function getBoundingSphere(grid: Grid, boundingSphere?: Sphere3D) {
        if (!boundingSphere) boundingSphere = Sphere3D();

        const dimensions = grid.cells.space.dimensions as Vec3;
        const transform = Grid.getGridToCartesianTransform(grid);
        return Sphere3D.fromDimensionsAndTransform(boundingSphere, dimensions, transform);
    }

    /**
     * Compute histogram with given bin count.
     * Cached on the Grid object.
     */
    export function getHistogram(grid: Grid, binCount: number) {
        let histograms = (grid as any)._historams as { [binCount: number]: Histogram };
        if (!histograms) {
            histograms = (grid as any)._historams = { };
        }
        if (!histograms[binCount]) {
            histograms[binCount] = calculateHistogram(grid.cells.data, binCount, { min: grid.stats.min, max: grid.stats.max });
        }
        return histograms[binCount];
    }

    export function makeGetTrilinearlyInterpolated(grid: Grid, transform: 'none' | 'relative') {
        const cartnToGrid = Grid.getGridToCartesianTransform(grid);
        Mat4.invert(cartnToGrid, cartnToGrid);
        const gridCoords = Vec3();

        const { stats } = grid;
        const { dimensions, get } = grid.cells.space;
        const data = grid.cells.data;

        const [mi, mj, mk] = dimensions;
        const getValue = (i: number, j: number, k: number) => get(data, i, j, k);

        return function getTrilinearlyInterpolated(position: Vec3): number {
            v3transformMat4(gridCoords, position, cartnToGrid);

            const value = trilinearlyInterpolate(gridCoords, mi, mj, mk, getValue);
            if (Number.isNaN(value)) return value;

            if (transform === 'relative') {
                return (value - stats.mean) / stats.sigma;
            } else {
                return value;
            }
        };
    }

    /**
     * Core trilinear interpolation function.
     * @param gridCoords - Position in grid coordinates (fractional indices)
     * @param mi, mj, mk - Grid dimensions
     * @param getValue - Function to get value at integer grid coordinates
     * @returns Interpolated value or NaN if out of bounds
     */
    export function trilinearlyInterpolate(
        gridCoords: Vec3,
        mi: number, mj: number, mk: number,
        getValue: (i: number, j: number, k: number) => number
    ): number {
        const i = Math.trunc(gridCoords[0]);
        const j = Math.trunc(gridCoords[1]);
        const k = Math.trunc(gridCoords[2]);

        if (i < 0 || i >= mi || j < 0 || j >= mj || k < 0 || k >= mk) {
            return Number.NaN;
        }

        const u = gridCoords[0] - i;
        const v = gridCoords[1] - j;
        const w = gridCoords[2] - k;

        const ii = Math.min(i + 1, mi - 1);
        const jj = Math.min(j + 1, mj - 1);
        const kk = Math.min(k + 1, mk - 1);

        let a = getValue(i, j, k);
        let b = getValue(ii, j, k);
        let c = getValue(i, jj, k);
        let d = getValue(ii, jj, k);
        const x = lerp(lerp(a, b, u), lerp(c, d, u), v);

        a = getValue(i, j, kk);
        b = getValue(ii, j, kk);
        c = getValue(i, jj, kk);
        d = getValue(ii, jj, kk);
        const y = lerp(lerp(a, b, u), lerp(c, d, u), v);

        return lerp(x, y, w);
    }

    const _a = Vec3(), _b = Vec3(), _c = Vec3(), _d = Vec3();
    const _ab = Vec3(), _cd = Vec3(), _x = Vec3(), _y = Vec3();

    /**
     * Core trilinear interpolation function for Vec3 values.
     * More efficient for interleaved data like gradients since getValue is called once per grid point.
     * @param gridCoords - Position in grid coordinates (fractional indices)
     * @param mi, mj, mk - Grid dimensions
     * @param getValue - Function to get Vec3 value at integer grid coordinates
     * @param out - Output Vec3
     * @returns true if interpolation succeeded, false if out of bounds
     */
    export function trilinearlyInterpolateVec3(
        gridCoords: Vec3,
        mi: number, mj: number, mk: number,
        getValue: (i: number, j: number, k: number, out: Vec3) => void,
        out: Vec3
    ): boolean {
        const i = Math.trunc(gridCoords[0]);
        const j = Math.trunc(gridCoords[1]);
        const k = Math.trunc(gridCoords[2]);

        if (i < 0 || i >= mi || j < 0 || j >= mj || k < 0 || k >= mk) {
            return false;
        }

        const u = gridCoords[0] - i;
        const v = gridCoords[1] - j;
        const w = gridCoords[2] - k;

        const ii = Math.min(i + 1, mi - 1);
        const jj = Math.min(j + 1, mj - 1);
        const kk = Math.min(k + 1, mk - 1);

        // Interpolate in the k plane
        getValue(i, j, k, _a);
        getValue(ii, j, k, _b);
        getValue(i, jj, k, _c);
        getValue(ii, jj, k, _d);
        v3lerp(_ab, _a, _b, u);
        v3lerp(_cd, _c, _d, u);
        v3lerp(_x, _ab, _cd, v);

        // Interpolate in the k+1 plane
        getValue(i, j, kk, _a);
        getValue(ii, j, kk, _b);
        getValue(i, jj, kk, _c);
        getValue(ii, jj, kk, _d);
        v3lerp(_ab, _a, _b, u);
        v3lerp(_cd, _c, _d, u);
        v3lerp(_y, _ab, _cd, v);

        // Final interpolation between planes
        v3lerp(out, _x, _y, w);
        return true;
    }

    type Gradients = {
        values: Float32Array,
        magnitude: {
            min: number,
            max: number,
            mean: number,
            sigma: number
        }
    };

    /**
     * Pre-compute gradients at each grid cell using central differences.
     * Returns a single Float32Array with interleaved xyz components (x1, y1, z1, x2, y2, z2, ...).
     * Cached on the Grid object.
     */
    export function getGradients(grid: Grid): Gradients {
        if ((grid as any)._gradients) {
            return (grid as any)._gradients as Gradients;
        }

        const gradients = computeGradients(grid);
        (grid as any)._gradients = gradients;
        return gradients;
    }

    function computeGradients(grid: Grid): Gradients {
        const { dimensions, get, dataOffset } = grid.cells.space;
        const data = grid.cells.data;
        const [mi, mj, mk] = dimensions;

        const n = mi * mj * mk;
        const values = new Float32Array(n * 3);

        let min = Infinity;
        let max = -Infinity;
        let sum = 0;
        let sumSq = 0;

        for (let k = 0; k < mk; ++k) {
            for (let j = 0; j < mj; ++j) {
                for (let i = 0; i < mi; ++i) {
                    const idx = dataOffset(i, j, k) * 3;

                    // Use central differences where possible, forward/backward at boundaries
                    const im = Math.max(0, i - 1);
                    const ip = Math.min(mi - 1, i + 1);
                    const jm = Math.max(0, j - 1);
                    const jp = Math.min(mj - 1, j + 1);
                    const km = Math.max(0, k - 1);
                    const kp = Math.min(mk - 1, k + 1);

                    // Gradient components (using central differences with proper divisor)
                    const gx = (get(data, ip, j, k) - get(data, im, j, k)) / (ip - im || 1);
                    const gy = (get(data, i, jp, k) - get(data, i, jm, k)) / (jp - jm || 1);
                    const gz = (get(data, i, j, kp) - get(data, i, j, km)) / (kp - km || 1);

                    values[idx] = gx;
                    values[idx + 1] = gy;
                    values[idx + 2] = gz;

                    const mag = gx * gx + gy * gy + gz * gz;
                    if (mag < min) min = mag;
                    if (mag > max) max = mag;
                    sum += mag;
                    sumSq += mag * mag;
                }
            }
        }

        if (min === Infinity) min = 0;
        if (max === -Infinity) max = 1;
        min = Math.sqrt(min);
        max = Math.sqrt(max);

        const mean = sum / n;
        const sigma = Math.sqrt(sumSq / n);

        return { values, magnitude: { min, max, mean, sigma } };
    }

    /**
     * Create a function that returns trilinearly interpolated gradient at a grid position.
     * The gradient is pre-computed at grid cells and interpolated for smooth results.
     */
    export function makeGetInterpolatedGradient(grid: Grid) {
        const { values: g } = getGradients(grid);
        const { dimensions, dataOffset } = grid.cells.space;
        const [mi, mj, mk] = dimensions;

        const getGradientVec3 = (i: number, j: number, k: number, out: Vec3) => {
            const idx = dataOffset(i, j, k) * 3;
            out[0] = g[idx];
            out[1] = g[idx + 1];
            out[2] = g[idx + 2];
        };

        return function getInterpolatedGradient(gridCoords: Vec3, out: Vec3): boolean {
            return trilinearlyInterpolateVec3(gridCoords, mi, mj, mk, getGradientVec3, out);
        };
    }
}

export { Grid };