/**
 * Copyright (c) 2018-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from '../../../../mol-model/structure';
import { Task } from '../../../../mol-task';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { GaussianDensityTexture, GaussianDensityTexture2d } from '../../../../mol-math/geometry/gaussian-density/gpu';
import { Texture } from '../../../../mol-gl/webgl/texture';
import { WebGLContext } from '../../../../mol-gl/webgl/context';
import { getUnitConformationAndRadius, getStructureConformationAndRadius, CommonSurfaceParams, ensureReasonableResolution } from './common';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';
import { GaussianDensityCPU } from '../../../../mol-math/geometry/gaussian-density/cpu';
import { SizeTheme } from '../../../../mol-theme/size';

export const GaussianDensityParams = {
    resolution: PD.Numeric(1, { min: 0.1, max: 20, step: 0.1 }, { description: 'Grid resolution/cell spacing.', ...BaseGeometry.CustomQualityParamInfo }),
    radiusOffset: PD.Numeric(0, { min: 0, max: 10, step: 0.1 }, { description: 'Extra/offset radius added to the atoms/coarse elements for gaussian calculation. Useful to create coarse, low resolution surfaces.' }),
    smoothness: PD.Numeric(1.5, { min: 0.5, max: 2.5, step: 0.1 }, { description: 'Smoothness of the gausian surface, lower is smoother.' }),
    ...CommonSurfaceParams
};
export const DefaultGaussianDensityProps = PD.getDefaultValues(GaussianDensityParams);
export type GaussianDensityProps = typeof DefaultGaussianDensityProps

//

export function getTextureMaxCells(webgl: WebGLContext, structure?: Structure) {
    const d = webgl.maxTextureSize / 3;
    return (d * d) / Math.max(1, (structure ? structure.units.length / 16 : 1));
}

//

export function computeUnitGaussianDensity(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: GaussianDensityProps) {
    const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props);
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensityCPU(ctx, position, boundary.box, radius, p);
    });
}

export function computeUnitGaussianDensityTexture(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props, getTextureMaxCells(webgl, structure));
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture(webgl, position, boundary.box, radius, p, texture);
    });
}

export function computeUnitGaussianDensityTexture2d(structure: Structure, unit: Unit, sizeTheme: SizeTheme<any>, powerOfTwo: boolean, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props, getTextureMaxCells(webgl, structure));
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture2d(webgl, position, boundary.box, radius, powerOfTwo, p, texture);
    });
}

//

export function computeStructureGaussianDensity(structure: Structure, sizeTheme: SizeTheme<any>, props: GaussianDensityProps) {
    const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props);
    return Task.create('Gaussian Density', async ctx => {
        return await GaussianDensityCPU(ctx, position, boundary.box, radius, p);
    });
}

export function computeStructureGaussianDensityTexture(structure: Structure, sizeTheme: SizeTheme<any>, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props);
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture(webgl, position, boundary.box, radius, p, texture);
    });
}

export function computeStructureGaussianDensityTexture2d(structure: Structure, sizeTheme: SizeTheme<any>, powerOfTwo: boolean, props: GaussianDensityProps, webgl: WebGLContext, texture?: Texture) {
    const { box } = structure.lookup3d.boundary;
    const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
    const p = ensureReasonableResolution(boundary.box, props);
    return Task.create('Gaussian Density', async ctx => {
        return GaussianDensityTexture2d(webgl, position, box, radius, powerOfTwo, p, texture);
    });
}
