/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import type { Unit, Structure } from '@molstar/model/model/structure';
import { Task } from '@molstar/core/task';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { GaussianDensityTexture, GaussianDensityTexture2d } from '@molstar/graphics/geo/gaussian-density/gpu';
import { isComputeGaussianDensitySupported } from '@molstar/graphics/gl/compute/gaussian-density';
import type { ComputeContext } from '@molstar/graphics/gl/compute';
import { Box3D } from '@molstar/core/math/geometry';
import { Vec3 } from '@molstar/core/math/linear-algebra';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import {
  getUnitConformationAndRadius,
  getStructureConformationAndRadius,
  CommonSurfaceParams,
  ensureReasonableResolution,
} from './common.js';
import { BaseGeometry } from '@molstar/graphics/geo/geometry/base';
import { GaussianDensityCPU } from '@molstar/core/math/geometry/gaussian-density/cpu';
import type { SizeTheme } from '@molstar/graphics/theme/size';

export const GaussianDensityParams = {
  resolution: PD.Numeric(
    1,
    { min: 0.1, max: 20, step: 0.1 },
    { description: 'Grid resolution/cell spacing.', ...BaseGeometry.CustomQualityParamInfo },
  ),
  radiusOffset: PD.Numeric(
    0,
    { min: 0, max: 10, step: 0.1 },
    {
      description:
        'Extra/offset radius added to the atoms/coarse elements for gaussian calculation. Useful to create coarse, low resolution surfaces.',
    },
  ),
  smoothness: PD.Numeric(
    1.5,
    { min: 1, max: 3, step: 0.1 },
    { description: 'Smoothness of the gausian surface, lower is smoother.' },
  ),
  floodfill: PD.Select('off', PD.arrayToOptions(['off', 'inside', 'outside']), {
    description: 'If and how to floodfill the gaussian surface. Note that this disables GPU support.',
  }),
  ...CommonSurfaceParams,
};
export const DefaultGaussianDensityProps = PD.getDefaultValues(GaussianDensityParams);
export type GaussianDensityProps = typeof DefaultGaussianDensityProps;

//

/**
 * Cells of a 3d density computed with kernels (WebGPU), by a device limit as the 2d atlas of WebGL
 * (`getTextureMaxCells`): the marching cubes keep a u32 count and offset per cell in storage buffers, so a cell
 * count beyond `maxStorageBufferSize / 4` fails to extract. That is 2²⁵ cells with the core limits (128 MiB storage
 * bindings, a 256 MiB rgba fp16 density), more with the adapter's (requested by default), where the GPU memory, not
 * reported by WebGPU, is the practical limit
 */
export function getComputeDensityMaxCells(compute: ComputeContext) {
  return Math.floor(compute.maxStorageBufferSize / 4);
}

/**
 * As `ensureReasonableResolution`, for the 3d density computed with kernels: limited by all cells of the grid, not by
 * the slices of a 2d atlas (sized by the 2d texture limit, which is smaller on WebGPU)
 */
function ensureComputeResolution<T>(box: Box3D, props: { resolution: number } & T, compute: ComputeContext) {
  const size = Box3D.size(Vec3(), box);
  const minResolution = Math.cbrt((size[0] * size[1] * size[2]) / getComputeDensityMaxCells(compute));
  return props.resolution < minResolution ? { ...props, resolution: minResolution } : props;
}

export function getTextureMaxCells(graphics: GraphicsContext, structure?: Structure) {
  const d = graphics.resources.limits.maxTextureSize / 3;
  return (d * d) / Math.max(1, structure ? structure.units.length / 16 : 1);
}

//

export function computeUnitGaussianDensity(
  structure: Structure,
  unit: Unit,
  sizeTheme: SizeTheme<any>,
  props: GaussianDensityProps,
) {
  const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
  const p = ensureReasonableResolution(boundary.box, props);
  return Task.create('Gaussian Density', async (ctx) => {
    return await GaussianDensityCPU(ctx, position, boundary.box, radius, p);
  });
}

export function computeUnitGaussianDensityTexture(
  structure: Structure,
  unit: Unit,
  sizeTheme: SizeTheme<any>,
  props: Omit<GaussianDensityProps, 'floodfill'>,
  graphics: GraphicsContext,
  texture?: Texture,
  growTexture?: boolean,
) {
  const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
  const p = isComputeGaussianDensitySupported(graphics)
    ? ensureComputeResolution(boundary.box, props, graphics.compute)
    : ensureReasonableResolution(boundary.box, props, getTextureMaxCells(graphics, structure));
  return GaussianDensityTexture(graphics, position, boundary.box, radius, p, texture, growTexture);
}

export function computeUnitGaussianDensityTexture2d(
  structure: Structure,
  unit: Unit,
  sizeTheme: SizeTheme<any>,
  powerOfTwo: boolean,
  props: GaussianDensityProps,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  const { position, boundary, radius } = getUnitConformationAndRadius(structure, unit, sizeTheme, props);
  const p = ensureReasonableResolution(boundary.box, props, getTextureMaxCells(graphics, structure));
  return GaussianDensityTexture2d(graphics, position, boundary.box, radius, powerOfTwo, p, texture);
}

//

export function computeStructureGaussianDensity(
  structure: Structure,
  sizeTheme: SizeTheme<any>,
  props: GaussianDensityProps,
) {
  const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
  const p = ensureReasonableResolution(boundary.box, props);
  return Task.create('Gaussian Density', async (ctx) => {
    return await GaussianDensityCPU(ctx, position, boundary.box, radius, p);
  });
}

export function computeStructureGaussianDensityTexture(
  structure: Structure,
  sizeTheme: SizeTheme<any>,
  props: Omit<GaussianDensityProps, 'floodfill'>,
  graphics: GraphicsContext,
  texture?: Texture,
  growTexture?: boolean,
) {
  const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
  const p = isComputeGaussianDensitySupported(graphics)
    ? ensureComputeResolution(boundary.box, props, graphics.compute)
    : ensureReasonableResolution(boundary.box, props);
  return GaussianDensityTexture(graphics, position, boundary.box, radius, p, texture, growTexture);
}

export function computeStructureGaussianDensityTexture2d(
  structure: Structure,
  sizeTheme: SizeTheme<any>,
  powerOfTwo: boolean,
  props: GaussianDensityProps,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  const { box } = structure.lookup3d.boundary;
  const { position, boundary, radius } = getStructureConformationAndRadius(structure, sizeTheme, props);
  const p = ensureReasonableResolution(boundary.box, props);
  return GaussianDensityTexture2d(graphics, position, box, radius, powerOfTwo, p, texture);
}
