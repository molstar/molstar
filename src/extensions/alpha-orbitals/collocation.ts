/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Inspired by https://github.com/dgasmith/gau2grid.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Vec3 } from '../../mol-math/linear-algebra';
import { RuntimeContext } from '../../mol-task';
import { arrayMin } from '../../mol-util/array';
import { AlphaOrbital, CubeGridInfo } from './data-model';
import { normalizeBasicOrder, SphericalFunctions } from './spherical-functions';

export async function sphericalCollocation(
    grid: CubeGridInfo,
    orbital: AlphaOrbital,
    taskCtx: RuntimeContext
) {
    const { basis, sphericalOrder, cutoffThreshold } = grid.params;
    let baseCount = 0;

    for (const atom of basis.atoms) {
        for (const shell of atom.shells) {
            for (const L of shell.angularMomentum) {
                if (L > 4) {
                    // TODO: will L > 4 be required? Would need to precompute more functions in that case.
                    throw new Error('Angular momentum L > 4 not supported.');
                }

                baseCount += 2 * L + 1;
            }
        }
    }

    const matrix = new Float32Array(grid.npoints);

    let baseIndex = 0;
    for (const atom of basis.atoms) {
        for (const shell of atom.shells) {
            let amIndex = 0;
            for (const L of shell.angularMomentum) {
                const alpha = normalizeBasicOrder(
                    L,
                    orbital.alpha.slice(baseIndex, baseIndex + 2 * L + 1),
                    sphericalOrder
                );
                baseIndex += 2 * L + 1;

                collocationBasis(
                    matrix,
                    grid,
                    L,
                    shell.coefficients[amIndex++],
                    shell.exponents,
                    atom.center as unknown as Vec3,
                    cutoffThreshold,
                    alpha
                );

                if (taskCtx.shouldUpdate) {
                    await taskCtx.update({
                        message: 'Computing...',
                        current: baseIndex,
                        max: baseCount,
                        isIndeterminate: false,
                    });
                }
            }
        }
    }

    return matrix;
}

function collocationBasis(
    matrix: Float32Array,
    grid: CubeGridInfo,
    L: number,
    coefficients: number[],
    exponents: number[],
    center: Vec3,
    cutoffThreshold: number,
    alpha: number[]
) {
    const ncoeff = exponents.length;
    const sphericalFunc = SphericalFunctions[L];

    const cx = center[0],
        cy = center[1],
        cz = center[2];
    const ny = grid.dimensions[1],
        nz = grid.dimensions[2];
    const gdx = grid.delta[0],
        gdy = grid.delta[1],
        gdz = grid.delta[2];
    const sx = grid.box.min[0],
        sy = grid.box.min[1],
        sz = grid.box.min[2];

    const cutoffRadius =
        cutoffThreshold > 0
            ? Math.sqrt(-Math.log(cutoffThreshold) / arrayMin(exponents))
            : 10000;
    const cutoffSquared = cutoffRadius * cutoffRadius;

    const radiusBox = getRadiusBox(grid, center, cutoffRadius);
    const iMin = radiusBox[0][0],
        jMin = radiusBox[0][1],
        kMin = radiusBox[0][2];
    const iMax = radiusBox[1][0],
        jMax = radiusBox[1][1],
        kMax = radiusBox[1][2];

    for (let i = iMin; i <= iMax; i++) {
        const x = sx + gdx * i - cx;
        const oX = i * ny * nz;

        for (let j = jMin; j <= jMax; j++) {
            const y = sy + gdy * j - cy;
            const oY = oX + j * nz;
            for (let k = kMin; k <= kMax; k++) {
                const z = sz + gdz * k - cz;
                const R2 = x * x + y * y + z * z;

                if (R2 > cutoffSquared) {
                    continue;
                }

                let gaussianSum = 0;
                for (let c = 0; c < ncoeff; c++) {
                    gaussianSum +=
                        coefficients[c] * Math.exp(-exponents[c] * R2);
                }

                const sphericalSum = L === 0 ? alpha[0] : sphericalFunc(alpha, x, y, z);


                matrix[k + oY] += gaussianSum * sphericalSum;
            }
        }
    }
}

function getRadiusBox(grid: CubeGridInfo, center: Vec3, radius: number) {
    const r = Vec3.create(radius, radius, radius);
    const min = Vec3.scaleAndAdd(Vec3(), center, r, -1);
    const max = Vec3.add(Vec3(), center, r);

    Vec3.sub(min, min, grid.box.min);
    Vec3.sub(max, max, grid.box.min);

    Vec3.div(min, min, grid.delta);
    Vec3.floor(min, min);
    Vec3.max(min, min, Vec3());

    Vec3.div(max, max, grid.delta);
    Vec3.ceil(max, max);
    Vec3.min(max, max, Vec3.subScalar(Vec3(), grid.dimensions, 1));

    return [min, max];
}
