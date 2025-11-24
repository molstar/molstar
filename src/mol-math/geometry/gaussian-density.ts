/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Box3D, DensityData } from '../geometry';
import { PositionData } from './common';
import { WebGLContext } from '../../mol-gl/webgl/context';
import { Texture } from '../../mol-gl/webgl/texture';
import { GaussianDensityTexture2d, GaussianDensityTexture3d } from './gaussian-density/gpu';
import { Task } from '../../mol-task/task';
import { GaussianDensityCPU } from './gaussian-density/cpu';
import { Mat4 } from '../linear-algebra/3d/mat4';
import { Vec3 } from '../linear-algebra/3d/vec3';
import { Vec2 } from '../linear-algebra/3d/vec2';

export const DefaultGaussianDensityProps = {
    resolution: 1,
    radiusOffset: 0,
    smoothness: 1.5,
};
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

export type GaussianDensityData = {
    radiusFactor: number
} & DensityData

export type GaussianDensityTextureData = {
    radiusFactor: number
    resolution: number
    maxRadius: number
    transform: Mat4,
    texture: Texture,
    bbox: Box3D,
    gridDim: Vec3,
    gridTexDim: Vec3
    gridDataDim: Vec3
    gridTexScale: Vec2
}

export function computeGaussianDensity(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps) {
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensityCPU(ctx, position, box, radius, props);
    });
}

export function computeGaussianDensityTexture(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    return _computeGaussianDensityTexture(webgl.isWebGL2 ? '3d' : '2d', position, box, radius, props, webgl, texture);
}

export function computeGaussianDensityTexture2d(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    return _computeGaussianDensityTexture('2d', position, box, radius, props, webgl, texture);
}

export function computeGaussianDensityTexture3d(position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    return _computeGaussianDensityTexture('2d', position, box, radius, props, webgl, texture);
}

function _computeGaussianDensityTexture(type: '2d' | '3d', position: PositionData, box: Box3D, radius: (index: number) => number, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    return Task.create('Gaussian Density', async ctx => {
        return type === '2d' ?
            GaussianDensityTexture2d(webgl, position, box, radius, false, props, texture) :
            GaussianDensityTexture3d(webgl, position, box, radius, props, texture);
    });
}