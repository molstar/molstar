/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '../geometry';
import { Vec3 } from '../linear-algebra';
import { RuntimeContext, Task } from 'mol-task';
import { PositionData, DensityData } from './common';
import { GaussianDensityCPU } from './gaussian-density/cpu';
import { WebGLContext } from 'mol-gl/webgl/context';

// import { GaussianDensityGPU } from './gaussian-density/gpu';
const GaussianDensityGPU = typeof document !== 'undefined'
    ? (require('./gaussian-density/gpu') as typeof import('./gaussian-density/gpu')).GaussianDensityGPU
    : void 0;

export const DefaultGaussianDensityProps = {
    resolution: 1,
    radiusOffset: 0,
    smoothness: 1.5,
    useGpu: true,
    webgl: undefined as WebGLContext | undefined
}
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

export function getDelta(box: Box3D, resolution: number) {
    const extent = Vec3.sub(Vec3.zero(), box.max, box.min)
    const size = Vec3.zero()
    Vec3.ceil(size, Vec3.scale(size, extent, resolution))
    const delta = Vec3.div(Vec3.zero(), extent, size)
    return delta
}

export function computeGaussianDensity(position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps) {
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, box, radius, props)
    });
}

export async function GaussianDensity(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps): Promise<DensityData> {
    if (props.useGpu) {
        if (!GaussianDensityGPU) throw 'GPU computation not supported on this platform';
        return await GaussianDensityGPU(ctx, position, box, radius, props)
    } else {
        return await GaussianDensityCPU(ctx, position, box, radius, props)
    }
}