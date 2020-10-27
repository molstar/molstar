/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Inspired by https://github.com/dgasmith/gau2grid.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArray } from '../../mol-data/util';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Box3D } from '../../mol-math/geometry';
import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { Grid } from '../../mol-model/volume';
import { getNonStandardResidueQueries } from '../../mol-plugin-state/helpers/structure-selection-query';
import { Task } from '../../mol-task';
import { arrayMax, arrayMin, arrayRms } from '../../mol-util/array';
import { CollocationParams, sphericalCollocation } from './collocation';
import { AlphaOrbitalsPass } from './gpu/pass';
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
    params: SphericalCollocationParams, webgl?: WebGLContext
): Task<CubeGrid> {
    return Task.create('Spherical Collocation Grid', async (ctx) => {
        const centers = params.basis.atoms.map(a => a.center);

        const cParams: CollocationParams = {
            grid: initBox(centers, params.gridSpacing, params.boxExpand),
            basis: params.basis,
            alphaOrbitals: params.alphaOrbitals,
            cutoffThreshold: params.cutoffThreshold,
            sphericalOrder: params.sphericalOrder
        };

        console.log(cParams);

        console.time('gpu');
        const pass = new AlphaOrbitalsPass(webgl!, cParams);
        const matrixGL = pass.getData();
        console.timeEnd('gpu');

        // TODO: remove the 2nd run
        console.time('gpu');
        const pass0 = new AlphaOrbitalsPass(webgl!, cParams);
        pass0.getData();
        console.timeEnd('gpu');

        // if (false && webgl) {
        // } else {
        // console.time('cpu');
        // const matrix = await sphericalCollocation(cParams, ctx);
        // console.timeEnd('cpu');
        // // }

        console.log(matrixGL);
        // console.log(matrix);

        // for (let i = 0; i < matrixGL.length; i++) {
        //     if (Math.abs(matrixGL[i] - matrix[i]) > 1e-4) {
        //         console.log('err', i, matrixGL[i], matrix[i]);
        //         // console.log()
        //         break;
        //     }
        // }

        return createCubeGrid(cParams.grid, matrixGL, [0, 1, 2]);
    });
}

const BohrToAngstromFactor = 0.529177210859;

function createCubeGrid(gridInfo: CubeGridInfo, values: Float32Array, axisOrder: number[]) {
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
            Tensor.Space(gridInfo.dimensions, axisOrder, Float32Array),
            (values as any) as Tensor.Data
        ),
        stats: {
            min: arrayMin(values),
            max: arrayMax(values),
            mean: arrayMax(values),
            sigma: arrayRms(values),
        },
    };

    // TODO: when using GPU rendering, the cumulative sum can be computed
    // along the ray on the fly
    console.time('iso');
    const isovalues = computeIsocontourValues(values, 0.85);
    console.timeEnd('iso');

    console.log(isovalues);

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
    input: Float32Array,
    cumulativeThreshold: number
) {
    let weightSum = 0;
    for (let i = 0, _i = input.length; i < _i; i++) {
        const v = input[i];
        const w = v * v;
        weightSum += w;
    }
    const avgWeight = weightSum / input.length;
    let minWeight = 3 * avgWeight;

    // do not try to identify isovalues for degenerate data
    // e.g. all values are almost same
    if (Math.abs(avgWeight - input[0] * input[0]) < 1e-5) {
        return { negative: void 0, positive: void 0 };
    }

    let size = 0;
    while (true) {
        let csum = 0;
        size = 0;
        for (let i = 0, _i = input.length; i < _i; i++) {
            const v = input[i];
            const w = v * v;
            if (w >= minWeight) {
                csum += w;
                size++;
            }
        }

        if (csum / weightSum > cumulativeThreshold) {
            break;
        }

        minWeight -= avgWeight;
    }

    const values = new Float32Array(size);
    const weights = new Float32Array(size);
    const indices = new Int32Array(size);

    let o = 0;
    for (let i = 0, _i = input.length; i < _i; i++) {
        const v = input[i];
        const w = v * v;
        if (w >= minWeight) {
            values[o] = v;
            weights[o] = w;
            indices[o] = o;
            o++;
        }
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