/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 *
 * Reading a per-voxel label array from isosurface vertex positions — shared by the color
 * themes of the volume tools.
 */

import { Mat4, Vec3 } from '../../mol-math/linear-algebra';
import { Grid, Volume } from '../../mol-model/volume';

function clampCell(v: number, n: number) {
    return Math.min(Math.max(v, 0), Math.max(n - 2, 0));
}

/**
 * Maps a world position to the label of the densest corner of the grid cell containing it.
 * Isosurface vertices lie on cell edges between an inside and an outside voxel, so taking the
 * densest corner always picks the inside voxel and surfaces are colored without speckle.
 */
export function makeLabelAtPosition(volume: Volume, labels: Uint8Array): (position: Vec3) => number {
    const { space, data } = volume.grid.cells;
    const values = data as unknown as ArrayLike<number>;
    const [nx, ny, nz] = space.dimensions as [number, number, number];
    const c2g = Mat4.invert(Mat4(), Grid.getGridToCartesianTransform(volume.grid));
    const g = Vec3();

    return (position: Vec3) => {
        Vec3.transformMat4(g, position, c2g);
        const i0 = clampCell(Math.floor(g[0]), nx), i1 = Math.min(i0 + 1, nx - 1);
        const j0 = clampCell(Math.floor(g[1]), ny), j1 = Math.min(j0 + 1, ny - 1);
        const k0 = clampCell(Math.floor(g[2]), nz), k1 = Math.min(k0 + 1, nz - 1);

        let best = -Infinity, label = 0;
        for (let i = i0; i <= i1; i++) {
            for (let j = j0; j <= j1; j++) {
                for (let k = k0; k <= k1; k++) {
                    const o = space.dataOffset(i, j, k);
                    const v = values[o];
                    if (v > best) {
                        best = v;
                        label = labels[o];
                    }
                }
            }
        }
        return label;
    };
}
