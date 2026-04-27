/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { Tensor } from '../../../mol-math/linear-algebra/tensor';

type Cells = { [i: number]: number; length: number };

/**
 * Mirror the volume along its first logical dimension (X / axis 0).
 * Swaps voxel pairs in-place: voxel[ix] ↔ voxel[nx-1-ix] for each (iy, iz) slice.
 */
export function flipVolumeX(cells: Cells, nx: number, ny: number, nz: number, space: Tensor.Space) {
    const half = Math.floor(nx / 2);
    for (let iz = 0; iz < nz; iz++) {
        for (let iy = 0; iy < ny; iy++) {
            for (let ix = 0; ix < half; ix++) {
                const a = space.dataOffset(ix, iy, iz);
                const b = space.dataOffset(nx - 1 - ix, iy, iz);
                const tmp = cells[a];
                cells[a] = cells[b];
                cells[b] = tmp;
            }
        }
    }
}

/**
 * Remove isolated blobs smaller than `minVoxels` by zeroing them out.
 *
 * Voxels with value >= `threshold` are treated as signal; connected components
 * are found via 6-connected iterative DFS. Components smaller than `minVoxels`
 * are zeroed in `cells`. Background (below threshold) voxels are left untouched.
 *
 * Returns the number of voxels zeroed.
 */
export function removeDust(cells: Cells, nx: number, ny: number, nz: number, space: Tensor.Space, minVoxels: number, threshold: number): number {
    const visited = new Uint8Array(cells.length);
    // Stack stores (ix, iy, iz) as consecutive triples — avoids needing getCoords.
    const stack: number[] = [];
    const component: number[] = [];
    let zeroed = 0;

    for (let iz = 0; iz < nz; iz++) {
        for (let iy = 0; iy < ny; iy++) {
            for (let ix = 0; ix < nx; ix++) {
                const off = space.dataOffset(ix, iy, iz);
                if (visited[off]) continue;
                visited[off] = 1;
                if (cells[off] < threshold) continue;

                // Found an unvisited signal voxel — flood-fill its component.
                component.length = 0;
                component.push(off);
                stack.length = 0;
                stack.push(ix, iy, iz);

                while (stack.length > 0) {
                    const cz = stack.pop()!;
                    const cy = stack.pop()!;
                    const cx = stack.pop()!;

                    const tryNeighbor = (nx2: number, ny2: number, nz2: number) => {
                        const noff = space.dataOffset(nx2, ny2, nz2);
                        if (visited[noff] || cells[noff] < threshold) return;
                        visited[noff] = 1;
                        component.push(noff);
                        stack.push(nx2, ny2, nz2);
                    };

                    if (cx > 0) tryNeighbor(cx - 1, cy, cz);
                    if (cx < nx - 1) tryNeighbor(cx + 1, cy, cz);
                    if (cy > 0) tryNeighbor(cx, cy - 1, cz);
                    if (cy < ny - 1) tryNeighbor(cx, cy + 1, cz);
                    if (cz > 0) tryNeighbor(cx, cy, cz - 1);
                    if (cz < nz - 1) tryNeighbor(cx, cy, cz + 1);
                }

                if (component.length < minVoxels) {
                    for (const o of component) cells[o] = 0;
                    zeroed += component.length;
                }
            }
        }
    }

    return zeroed;
}
