/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Mat4, Tensor, Vec3 } from '../../mol-math/linear-algebra';
import { Grid } from '../../mol-model/volume';
import { SphericalBasisOrder } from './spherical-functions';
import { Box3D, RegularGrid3d } from '../../mol-math/geometry';
import { arrayMin, arrayMax, arrayRms, arrayMean } from '../../mol-util/array';
import { ModelFormat } from '../../mol-model-formats/format';

// Note: generally contracted gaussians are currently not supported.
export interface SphericalElectronShell {
    exponents: number[];
    angularMomentum: number[];
    // number[] for each angular momentum
    coefficients: number[][];
}

export interface Basis {
    atoms: {
        // in Bohr units!
        center: [number, number, number];
        shells: SphericalElectronShell[];
    }[];
}

export interface AlphaOrbital {
    energy: number;
    occupancy: number;
    alpha: number[];
}

export interface CubeGridComputationParams {
    basis: Basis;
    /**
     * for each electron shell compute a cutoff radius as
     *    const cutoffRadius = Math.sqrt(-Math.log(cutoffThreshold) / arrayMin(exponents));
     */
    cutoffThreshold: number;
    sphericalOrder: SphericalBasisOrder;
    boxExpand: number;
    gridSpacing: number | [atomCountThreshold: number, spacing: number][];
    doNotComputeIsovalues?: boolean;
}

export interface CubeGridInfo {
    params: CubeGridComputationParams;
    dimensions: Vec3;
    box: Box3D;
    size: Vec3;
    npoints: number;
    delta: Vec3;
}

export interface CubeGrid {
    grid: Grid;
    isovalues?: { negative?: number; positive?: number };
}

export type CubeGridFormat = ModelFormat<CubeGrid>;

export function CubeGridFormat(grid: CubeGrid): CubeGridFormat {
    return { name: 'custom grid', kind: 'cube-grid', data: grid };
}

export function isCubeGridData(f: ModelFormat): f is CubeGridFormat {
    return f.kind === 'cube-grid';
}

export function initCubeGrid(params: CubeGridComputationParams): CubeGridInfo {
    const geometry = params.basis.atoms.map(a => a.center);
    const { gridSpacing: spacing, boxExpand: expand } = params;

    const count = geometry.length;
    const box = Box3D.expand(
        Box3D(),
        Box3D.fromVec3Array(Box3D(), geometry as unknown as Vec3[]),
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
        params,
        box,
        dimensions,
        size,
        npoints: dimensions[0] * dimensions[1] * dimensions[2],
        delta: Vec3.div(Vec3(), size, Vec3.subScalar(Vec3(), dimensions, 1)),
    };
}

const BohrToAngstromFactor = 0.529177210859;

export function createGrid(gridInfo: RegularGrid3d, values: Float32Array, axisOrder: number[]) {
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
            mean: arrayMean(values),
            sigma: arrayRms(values),
        },
    };

    return grid;
}
