/**
 * Copyright (c) 2018-2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../../../mol-model/structure';
import { GaussianDensity } from '../../../../mol-math/geometry/gaussian-density';
import { Task } from '../../../../mol-task';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { GaussianDensityTexture, GaussianDensityTexture2d } from '../../../../mol-math/geometry/gaussian-density/gpu';
import { Texture } from '../../../../mol-gl/webgl/texture';
import { WebGLContext } from '../../../../mol-gl/webgl/context';
import { getUnitConformationAndRadius, getStructureConformationAndRadius } from './common';

export const GaussianDensityParams = {
    resolution: PD.Numeric(1, { min: 0.1, max: 20, step: 0.1 }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    smoothness: PD.Numeric(1.5, { min: 0.5, max: 2.5, step: 0.1 }),
    useGpu: PD.Boolean(false),
    ignoreHydrogens: PD.Boolean(false),
}
export const DefaultGaussianDensityProps = PD.getDefaultValues(GaussianDensityParams)
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

export const GaussianDensityTextureParams = {
    resolution: PD.Numeric(1, { min: 0.1, max: 20, step: 0.1 }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }),
    smoothness: PD.Numeric(1.5, { min: 0.5, max: 2.5, step: 0.1 }),
    ignoreHydrogens: PD.Boolean(false),
}
export const DefaultGaussianDensityTextureProps = PD.getDefaultValues(GaussianDensityTextureParams)
export type GaussianDensityTextureProps = typeof DefaultGaussianDensityTextureProps

//

export function computeUnitGaussianDensity(unit: Unit, props: GaussianDensityProps, webgl?: WebGLContext) {
    const { position, radius } = getUnitConformationAndRadius(unit, props.ignoreHydrogens)
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, unit.lookup3d.boundary.box, radius, props, webgl);
    });
}

export function computeUnitGaussianDensityTexture(unit: Unit, props: GaussianDensityTextureProps, webgl: WebGLContext, texture?: Texture) {
    const { position, radius } = getUnitConformationAndRadius(unit, props.ignoreHydrogens)
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture(webgl, position, unit.lookup3d.boundary.box, radius, props, texture);
    });
}

export function computeUnitGaussianDensityTexture2d(unit: Unit, props: GaussianDensityTextureProps, webgl: WebGLContext, texture?: Texture) {
    const { position, radius } = getUnitConformationAndRadius(unit, props.ignoreHydrogens)
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture2d(webgl, position, unit.lookup3d.boundary.box, radius, props, texture);
    });
}

//

export function computeStructureGaussianDensity(structure: Structure, props: GaussianDensityProps, webgl?: WebGLContext) {
    const { position, radius } = getStructureConformationAndRadius(structure, props.ignoreHydrogens)
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, structure.lookup3d.boundary.box, radius, props, webgl);
    });
}

export function computeStructureGaussianDensityTexture(structure: Structure, props: GaussianDensityTextureProps, webgl: WebGLContext, texture?: Texture) {
    const { position, radius } = getStructureConformationAndRadius(structure, props.ignoreHydrogens)
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture(webgl, position, structure.lookup3d.boundary.box, radius, props, texture);
    });
}