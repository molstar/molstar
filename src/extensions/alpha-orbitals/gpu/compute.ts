/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createGrid3dComputeRenderable } from '../../../mol-gl/compute/grid3d';
import { TextureSpec, UniformSpec } from '../../../mol-gl/renderable/schema';
import { WebGLContext } from '../../../mol-gl/webgl/context';
import { RuntimeContext } from '../../../mol-task';
import { ValueCell } from '../../../mol-util';
import { AlphaOrbital, CubeGridInfo } from '../data-model';
import { createTextureData, getNormalizedAlpha } from './data';
import { MAIN, UTILS } from './shader.frag';

const Schema = {
    tCenters: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    tInfo: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    tCoeff: TextureSpec('image-float32', 'rgb', 'float', 'nearest'),
    tAlpha: TextureSpec('image-float32', 'alpha', 'float', 'nearest'),
    uNCenters: UniformSpec('i'),
    uNAlpha: UniformSpec('i'),
    uNCoeff: UniformSpec('i'),
    uMaxCoeffs: UniformSpec('i'),
};

const Orbitals = createGrid3dComputeRenderable({
    schema: Schema,
    loopBounds: ['uNCenters', 'uMaxCoeffs'],
    mainCode: MAIN,
    utilCode: UTILS,
    returnCode: 'v',
    values(params: { grid: CubeGridInfo, orbital: AlphaOrbital }) {
        return createTextureData(params.grid, params.orbital);
    }
});

const Density = createGrid3dComputeRenderable({
    schema: {
        ...Schema,
        uOccupancy: UniformSpec('f'),
    },
    loopBounds: ['uNCenters', 'uMaxCoeffs'],
    mainCode: MAIN,
    utilCode: UTILS,
    returnCode: 'current + uOccupancy * v * v',
    values(params: { grid: CubeGridInfo, orbitals: AlphaOrbital[] }) {
        return {
            ...createTextureData(params.grid, params.orbitals[0]),
            uOccupancy: 0
        };
    },
    cumulative: {
        states(params: { grid: CubeGridInfo, orbitals: AlphaOrbital[] }) {
            return params.orbitals.filter(o => o.occupancy !== 0);
        },
        update({ grid }, state: AlphaOrbital, values) {
            const alpha = getNormalizedAlpha(grid.params.basis, state.alpha, grid.params.sphericalOrder);
            ValueCell.updateIfChanged(values.uOccupancy, state.occupancy);
            ValueCell.update(values.tAlpha, { width: alpha.length, height: 1, array: alpha });
        }
    }
});

export function gpuComputeAlphaOrbitalsGridValues(ctx: RuntimeContext, webgl: WebGLContext, grid: CubeGridInfo, orbital: AlphaOrbital) {
    return Orbitals(ctx, webgl, grid, { grid, orbital });
}

export function gpuComputeAlphaOrbitalsDensityGridValues(ctx: RuntimeContext, webgl: WebGLContext, grid: CubeGridInfo, orbitals: AlphaOrbital[]) {
    return Density(ctx, webgl, grid, { grid, orbitals });
}
