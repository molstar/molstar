/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D } from '../geometry';
import { RuntimeContext, Task } from '../../mol-task';
import { PositionData, DensityData } from './common';
import { GaussianDensityCPU } from './gaussian-density/cpu';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { GaussianDensityTexture2d, GaussianDensityTexture3d } from './gaussian-density/gpu';

// import { GaussianDensityGPU, GaussianDensityTexture } from './gaussian-density/gpu';
const GaussianDensityGPU = typeof document !== 'undefined'
    ? (require('./gaussian-density/gpu') as typeof import('./gaussian-density/gpu')).GaussianDensityGPU
    : void 0;
const GaussianDensityTexture = typeof document !== 'undefined'
    ? (require('./gaussian-density/gpu') as typeof import('./gaussian-density/gpu')).GaussianDensityTexture
    : void 0;

export const DefaultGaussianDensityGPUProps = {
    resolution: 1,
    radiusOffset: 0,
    smoothness: 1.5,
};
export type GaussianDensityGPUProps = typeof DefaultGaussianDensityGPUProps

export const DefaultGaussianDensityProps = {
    ...DefaultGaussianDensityGPUProps,
    useGpu: true,
};
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

export function computeGaussianDensity(position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps, webgl?: WebGLContext) {
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, box, radius, props, webgl);
    });
}

export async function GaussianDensity(ctx: RuntimeContext, position: PositionData, box: Box3D, radius: (index: number) => number,  props: GaussianDensityProps, webgl?: WebGLContext): Promise<DensityData> {
    if (props.useGpu) {
        if (!GaussianDensityGPU) throw 'GPU computation not supported on this platform';
        if (!webgl) throw 'No WebGL context provided';
        return GaussianDensityGPU(position, box, radius, props, webgl);
    } else {
        return await GaussianDensityCPU(ctx, position, box, radius, props);
    }
}

export function computeGaussianDensityTexture(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, webgl: WebGLContext, texture?: Texture) {
    return _computeGaussianDensityTexture(webgl.isWebGL2 ? '3d' : '2d', position, box, radius, props, webgl, texture);
}

export function computeGaussianDensityTexture2d(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, webgl: WebGLContext, texture?: Texture) {
    return _computeGaussianDensityTexture('2d', position, box, radius, props, webgl, texture);
}

export function computeGaussianDensityTexture3d(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, webgl: WebGLContext, texture?: Texture) {
    return _computeGaussianDensityTexture('2d', position, box, radius, props, webgl, texture);
}

function _computeGaussianDensityTexture(type: '2d' | '3d', position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityGPUProps, webgl: WebGLContext, texture?: Texture) {
    if (!GaussianDensityTexture) throw 'GPU computation not supported on this platform';
    return Task.create('Gaussian Density', async ctx => {
        return type === '2d' ?
            GaussianDensityTexture2d(webgl, position, box, radius, props, texture) :
            GaussianDensityTexture3d(webgl, position, box, radius, props, texture);
    });
}