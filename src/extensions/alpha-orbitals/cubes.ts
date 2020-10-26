/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Inspired by https://github.com/dgasmith/gau2grid.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArray } from '../../mol-data/util';
import { Box3D } from '../../mol-math/geometry';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { Grid } from '../../mol-model/volume';
import { Task } from '../../mol-task';
import { arrayMax, arrayMin, arrayRms } from '../../mol-util/array';
import { sphericalCollocation } from './collocation';
import { SphericalBasisOrder } from './orbitals';

export interface CubeGridInfo {
    dimensions: Vec3;
    box: Box3D;
    size: Vec3;
    npoints: number;
    delta: Vec3;
}

export interface CubeGrid {
    grid: Grid;
    isovalues: { negative?: number; positive?: number };
}

export interface Basis {
    atoms: {
        // in Bohr units!
        center: Vec3;
        shells: SphericalElectronShell[];
    }[];
}

// Note: generally contracted gaussians are currently not supported.
export interface SphericalElectronShell {
    exponents: number[];
    angularMomentum: number[];
    // number[] for each angular momentum
    coefficients: number[][];
}

export interface SphericalCollocationParams {
    basis: Basis;
    /**
     * for each electron shell compute a cutoff radius as
     *
     *    const cutoffRadius = Math.sqrt(-Math.log(cutoffThreshold) / arrayMin(exponents));
     *
     */
    cutoffThreshold: number;
    sphericalOrder: SphericalBasisOrder;
    boxExpand: number;
    gridSpacing: number | [atomCountThreshold: number, spacing: number][];
    alphaOrbitals: number[];
}

export function createSphericalCollocationGrid(
    params: SphericalCollocationParams
): Task<CubeGrid> {
    return Task.create('Spherical Collocation Grid', async (ctx) => {
        const grid = initBox(
            params.basis.atoms.map((a) => a.center),
            params.gridSpacing,
            params.boxExpand
        );

        const matrix = await sphericalCollocation(
            {
                grid,
                basis: params.basis,
                alphaOrbitals: params.alphaOrbitals,
                cutoffThreshold: params.cutoffThreshold,
                sphericalOrder: 'cca-reverse', // renamed entos ordering
            },
            ctx
        );

        return createCubeGrid(grid, matrix);
    });
}

const BohrToAngstromFactor = 0.529177210859;

function createCubeGrid(gridInfo: CubeGridInfo, values: Float32Array) {
    const boxSize = Box3D.size(Vec3(), gridInfo.box);
    const boxOrigin = Vec3.clone(gridInfo.box.min);

    Vec3.scale(boxSize, boxSize, BohrToAngstromFactor);
    Vec3.scale(boxOrigin, boxOrigin, BohrToAngstromFactor);

    const scale = Mat4.fromScaling(
        Mat4(),
        Vec3.div(
            Vec3(),
            boxSize,
            Vec3.sub(Vec3(), gridInfo.dimensions, Vec3.create(1, 1, 1))
        )
    );
    const translate = Mat4.fromTranslation(Mat4(), boxOrigin);
    const matrix = Mat4.mul(Mat4(), translate, scale);

    const grid: Grid = {
        transform: { kind: 'matrix', matrix },
        cells: Tensor.create(
            Tensor.Space(gridInfo.dimensions, [0, 1, 2], Float32Array),
            (values as any) as Tensor.Data
        ),
        stats: {
            min: arrayMin(values),
            max: arrayMax(values),
            mean: arrayMax(values),
            sigma: arrayRms(values),
        },
    };

    const isovalues = computeIsocontourValues(values, 0.85);

    return { grid, isovalues };
}

function initBox(
    geometry: Vec3[],
    spacing: SphericalCollocationParams['gridSpacing'],
    expand: number
): CubeGridInfo {
    const count = geometry.length;
    const box = Box3D.expand(
        Box3D(),
        Box3D.fromVec3Array(Box3D(), geometry),
        Vec3.create(expand, expand, expand)
    );
    const size = Box3D.size(Vec3(), box);

    const spacingThresholds =
        typeof spacing === 'number' ? [[0, spacing]] : [...spacing];
    spacingThresholds.sort((a, b) => b[0] - a[0]);

    let s = 0.4;
    for (let i = 0; i <= spacingThresholds.length; i++) {
        s = spacingThresholds[i][1];
        if (spacingThresholds[i][0] <= count) break;
    }

    const dimensions = Vec3.ceil(Vec3(), Vec3.scale(Vec3(), size, 1 / s));
    return {
        box,
        dimensions,
        size,
        npoints: dimensions[0] * dimensions[1] * dimensions[2],
        delta: Vec3.div(Vec3(), size, Vec3.subScalar(Vec3(), dimensions, 1)),
    };
}

export function computeIsocontourValues(
    values: Float32Array,
    cumulativeThreshold: number
) {
    const size = values.length;
    const weights = new Float32Array(size);
    const indices = new Int32Array(size);

    let weightSum = 0;
    for (let i = 0; i < size; i++) {
        const v = values[i];
        const w = v * v;
        weights[i] = w;
        indices[i] = i;
        weightSum += w;
    }

    sortArray(
        indices,
        (indices, i, j) => weights[indices[j]] - weights[indices[i]]
    );

    let cweight = 0,
        cutoffIndex = 0;
    for (let i = 0; i < size; i++) {
        cweight += weights[indices[i]];

        if (cweight / weightSum >= cumulativeThreshold) {
            cutoffIndex = i;
            break;
        }
    }

    let positive = Number.POSITIVE_INFINITY,
        negative = Number.NEGATIVE_INFINITY;

    for (let i = 0; i < cutoffIndex; i++) {
        const v = values[indices[i]];
        if (v > 0) {
            if (v < positive) positive = v;
        } else if (v < 0) {
            if (v > negative) negative = v;
        }
    }

    return {
        negative: negative !== Number.NEGATIVE_INFINITY ? negative : void 0,
        positive: positive !== Number.POSITIVE_INFINITY ? positive : void 0,
    };
}
