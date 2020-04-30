/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { getUnitConformationAndRadius, getStructureConformationAndRadius, CommonSurfaceParams, ensureReasonableResolution } from './common';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';

const SharedGaussianDensityParams = {
    resolution: PD.Numeric(1, { min: 0.1, max: 20, step: 0.1 }, { description: 'Grid resolution/cell spacing.', ...BaseGeometry.CustomQualityParamInfo }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, { description: 'Extra/offset radius added to the atoms/coarse elements for gaussian calculation. Useful to create coarse, low resolution surfaces.' }),
    smoothness: PD.Numeric(1.5, { min: 0.5, max: 2.5, step: 0.1 }, { description: 'Smoothness of the gausian surface, lower is smoother.' }),
    ...CommonSurfaceParams
};

export const GaussianDensityParams = {
    ...SharedGaussianDensityParams,
    useGpu: PD.Boolean(false),
};
export const DefaultGaussianDensityProps = PD.getDefaultValues(GaussianDensityParams);
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

export const GaussianDensityTextureParams = {
    ...SharedGaussianDensityParams
};
export const DefaultGaussianDensityTextureProps = PD.getDefaultValues(GaussianDensityTextureParams);
export type GaussianDensityTextureProps = typeof DefaultGaussianDensityTextureProps

//

export function computeUnitGaussianDensity(structure: Structure, unit: Unit, props: GaussianDensityProps, webgl?: WebGLContext) {
    const { box } = unit.lookup3d.boundary;
    const p = ensureReasonableResolution(box, props);
    const { position, radius } = getUnitConformationAndRadius(structure, unit, p);
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, box, radius, p, webgl);
    });
}

export function computeUnitGaussianDensityTexture(structure: Structure, unit: Unit, props: GaussianDensityTextureProps, webgl: WebGLContext, texture?: Texture) {
    const { box } = unit.lookup3d.boundary;
    const p = ensureReasonableResolution(box, props);
    const { position, radius } = getUnitConformationAndRadius(structure, unit, p);
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture(webgl, position, box, radius, p, texture);
    });
}

export function computeUnitGaussianDensityTexture2d(structure: Structure, unit: Unit, props: GaussianDensityTextureProps, webgl: WebGLContext, texture?: Texture) {
    const { box } = unit.lookup3d.boundary;
    const p = ensureReasonableResolution(box, props);
    const { position, radius } = getUnitConformationAndRadius(structure, unit, p);
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture2d(webgl, position, box, radius, p, texture);
    });
}

//

export function computeStructureGaussianDensity(structure: Structure, props: GaussianDensityProps, webgl?: WebGLContext) {
    const { box } = structure.lookup3d.boundary;
    const p = ensureReasonableResolution(box, props);
    const { position, radius } = getStructureConformationAndRadius(structure, props.ignoreHydrogens, props.traceOnly);
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensity(ctx, position, box, radius, p, webgl);
    });
}

export function computeStructureGaussianDensityTexture(structure: Structure, props: GaussianDensityTextureProps, webgl: WebGLContext, texture?: Texture) {
    const { box } = structure.lookup3d.boundary;
    const p = ensureReasonableResolution(box, props);
    const { position, radius } = getStructureConformationAndRadius(structure, props.ignoreHydrogens, props.traceOnly);
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture(webgl, position, box, radius, p, texture);
    });
}