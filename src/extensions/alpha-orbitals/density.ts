/**
 * Copyright (c) 2020-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { sortArray } from '../../mol-data/util';
import { canComputeGrid3dOnGPU } from '../../mol-gl/compute/grid3d';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Task } from '../../mol-task';
import { isTimingMode } from '../../mol-util/debug';
import { AlphaOrbital, createGrid, CubeGrid, CubeGridComputationParams, initCubeGrid } from './data-model';
import { gpuComputeAlphaOrbitalsDensityGridValues } from './gpu/compute';

export function createSphericalCollocationDensityGrid(
    params: CubeGridComputationParams, orbitals: AlphaOrbital[], webgl?: WebGLContext
): Task<CubeGrid> {
    return Task.create('Spherical Collocation Grid', async (ctx) => {
        const cubeGrid = initCubeGrid(params);

        let matrix: Float32Array;
        if (canComputeGrid3dOnGPU(webgl)) {
            if (isTimingMode) webgl.timer.mark('createSphericalCollocationDensityGrid');
            matrix = await gpuComputeAlphaOrbitalsDensityGridValues(ctx, webgl, cubeGrid, orbitals);
            if (isTimingMode) webgl.timer.markEnd('createSphericalCollocationDensityGrid');
        } else {
            throw new Error('Missing OES_texture_float WebGL extension.');
        }

        const grid = createGrid(cubeGrid, matrix, [0, 1, 2]);
        let isovalues: { negative?: number, positive?: number } | undefined;

        if (!params.doNotComputeIsovalues) {
            isovalues = computeDensityIsocontourValues(matrix, 0.85);
        }

        return { grid, isovalues };
    });
}

export function computeDensityIsocontourValues(input: Float32Array, cumulativeThreshold: number) {
    let weightSum = 0;
    for (let i = 0, _i = input.length; i < _i; i++) {
        const v = input[i];
        const w = Math.abs(v);
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
            const w = Math.abs(v);
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
        const w = Math.abs(v);
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