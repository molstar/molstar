/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SpacegroupCell, Box3D } from '../../mol-math/geometry';
import { Tensor, Mat4, Vec3 } from '../../mol-math/linear-algebra';

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
}

export { Grid };