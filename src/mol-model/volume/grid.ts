/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SpacegroupCell, Box3D, Sphere3D } from '../../mol-math/geometry';
import { Tensor, Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Histogram, calculateHistogram } from '../../mol-math/histogram';
import { lerp } from '../../mol-math/interpolate';

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

    const _scale = Mat4.zero(), _translate = Mat4.zero();
    export function getGridToCartesianTransform(grid: Grid) {
        if (grid.transform.kind === 'matrix') {
            return Mat4.copy(Mat4(), grid.transform.matrix);
        }

        if (grid.transform.kind === 'spacegroup') {
            const { cells: { space } } = grid;
            const scale = Mat4.fromScaling(_scale, Vec3.div(Vec3.zero(), Box3D.size(Vec3.zero(), grid.transform.fractionalBox), Vec3.ofArray(space.dimensions)));
            const translate = Mat4.fromTranslation(_translate, grid.transform.fractionalBox.min);
            return Mat4.mul3(Mat4.zero(), grid.transform.cell.fromFractional, translate, scale);
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

        return function getTrilinearlyInterpolated(position: Vec3): number {
            Vec3.copy(gridCoords, position);
            Vec3.transformMat4(gridCoords, gridCoords, cartnToGrid);

            const i = Math.trunc(gridCoords[0]);
            const j = Math.trunc(gridCoords[1]);
            const k = Math.trunc(gridCoords[2]);

            if (i < 0 || i >= mi || j < 0 || j >= mj || k < 0 || k >= mk) {
                return Number.NaN;
            }

            const u = gridCoords[0] - i;
            const v = gridCoords[1] - j;
            const w = gridCoords[2] - k;

            // Tri-linear interpolation for the value
            const ii = Math.min(i + 1, mi - 1);
            const jj = Math.min(j + 1, mj - 1);
            const kk = Math.min(k + 1, mk - 1);

            let a = get(data, i, j, k);
            let b = get(data, ii, j, k);
            let c = get(data, i, jj, k);
            let d = get(data, ii, jj, k);
            const x = lerp(lerp(a, b, u), lerp(c, d, u), v);

            a = get(data, i, j, kk);
            b = get(data, ii, j, kk);
            c = get(data, i, jj, kk);
            d = get(data, ii, jj, kk);
            const y = lerp(lerp(a, b, u), lerp(c, d, u), v);

            const value = lerp(x, y, w);
            if (transform === 'relative') {
                return (value - stats.mean) / stats.sigma;
            } else {
                return value;
            }
        };
    }
}

export { Grid };