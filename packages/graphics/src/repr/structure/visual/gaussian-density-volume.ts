/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { VisualContext } from '../../visual.js';
import type { Structure, Unit } from '@molstar/model/model/structure';
import type { Theme } from '@molstar/graphics/theme/theme';
import {
  type GaussianDensityProps,
  computeStructureGaussianDensityTexture,
  computeUnitGaussianDensityTexture,
  GaussianDensityParams,
} from './util/gaussian.js';
import { DirectVolume } from '@molstar/graphics/geo/geometry/direct-volume/direct-volume';
import { ComplexDirectVolumeParams, type ComplexVisual, ComplexDirectVolumeVisual } from '../complex-visual.js';
import type { VisualUpdateState } from '../../util.js';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import {
  eachElement,
  eachSerialElement,
  ElementIterator,
  getElementLoci,
  getSerialElementLoci,
} from './util/element.js';
import { Sphere3D } from '@molstar/core/math/geometry';
import { UnitsDirectVolumeParams, type UnitsVisual, UnitsDirectVolumeVisual } from '../units-visual.js';
import { omitObjectKeys } from '@molstar/core/util/object';
import { isGaussianDensityTextureSupported } from '@molstar/graphics/geo/gaussian-density/gpu';

function createGaussianDensityVolume(
  ctx: VisualContext,
  structure: Structure,
  theme: Theme,
  props: Omit<GaussianDensityProps, 'floodfill'>,
  directVolume?: DirectVolume,
): DirectVolume {
  const graphics = ctx.graphics;
  if (!graphics || !isGaussianDensityTextureSupported(graphics)) {
    throw new Error('GaussianDensityVolume requires the GPU gaussian density');
  }

  const axisOrder = Vec3.create(0, 1, 2);
  const stats = { min: 0, max: 1, mean: 0.04, sigma: 0.01 };

  const create = (directVolume?: DirectVolume) => {
    // not one of the previous context (a backend switch or restore), which destroys it
    const previous = directVolume?.gridTexture.ref.value;
    const oldTexture = previous && graphics.resources.ownsTexture(previous) ? previous : undefined;
    const densityTextureData = computeStructureGaussianDensityTexture(
      structure,
      theme.size,
      props,
      graphics,
      oldTexture,
    );
    const { transform, texture, bbox, gridDim } = densityTextureData;

    const unitToCartn = Mat4.mul(Mat4(), transform, Mat4.fromScaling(Mat4(), gridDim));
    const cellDim = Mat4.getScaling(Vec3(), transform);

    const vol = DirectVolume.create(
      bbox,
      gridDim,
      transform,
      unitToCartn,
      cellDim,
      texture,
      stats,
      true,
      axisOrder,
      'byte',
      directVolume,
    );

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, densityTextureData.maxRadius);
    vol.setBoundingSphere(sphere);
    return vol;
  };

  return create(directVolume);
}

/**
 * Empty space skipping of the Gaussian volume: a sample this far ahead (in Å) is probed where the density is below 0.01
 * and the space between skipped if it is too (approximate: may step over thin features); volumes skip bricks the
 * transfer function makes transparent instead, which the GPU-made density has no CPU grid for
 */
const JumpLengthParam = PD.Numeric(4, { min: 0, max: 20, step: 0.1 });

export const GaussianDensityVolumeParams = {
  ...ComplexDirectVolumeParams,
  jumpLength: JumpLengthParam,
  ...omitObjectKeys(GaussianDensityParams, ['floodfill']),
  ignoreHydrogens: PD.Boolean(false),
  ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
  includeParent: PD.Boolean(false, { isHidden: true }),
};
export type GaussianDensityVolumeParams = typeof GaussianDensityVolumeParams;

export function GaussianDensityVolumeVisual(materialId: number): ComplexVisual<GaussianDensityVolumeParams> {
  return ComplexDirectVolumeVisual<GaussianDensityVolumeParams>(
    {
      defaultProps: PD.getDefaultValues(GaussianDensityVolumeParams),
      createGeometry: createGaussianDensityVolume,
      // the density packs the group with the value
      signatureDefines: () => ({ dPackedGroup: [true] }),
      createLocationIterator: ElementIterator.fromStructure,
      getLoci: getSerialElementLoci,
      eachLocation: eachSerialElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<GaussianDensityVolumeParams>,
        currentProps: PD.Values<GaussianDensityVolumeParams>,
      ) => {
        if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
        if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
        if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
        if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
        if (newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant) state.createGeometry = true;
        if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
        if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
      },
      dispose: (geometry: DirectVolume) => {
        geometry.gridTexture.ref.value.destroy();
      },
    },
    materialId,
  );
}

//

function createUnitsGaussianDensityVolume(
  ctx: VisualContext,
  unit: Unit,
  structure: Structure,
  theme: Theme,
  props: Omit<GaussianDensityProps, 'floodfill'>,
  directVolume?: DirectVolume,
): DirectVolume {
  const graphics = ctx.graphics;
  if (!graphics || !isGaussianDensityTextureSupported(graphics)) {
    throw new Error('GaussianDensityVolume requires the GPU gaussian density');
  }

  const axisOrder = Vec3.create(0, 1, 2);
  const stats = { min: 0, max: 1, mean: 0.04, sigma: 0.01 };

  const create = (directVolume?: DirectVolume) => {
    // not one of the previous context (a backend switch or restore), which destroys it
    const previous = directVolume?.gridTexture.ref.value;
    const oldTexture = previous && graphics.resources.ownsTexture(previous) ? previous : undefined;
    const densityTextureData = computeUnitGaussianDensityTexture(
      structure,
      unit,
      theme.size,
      props,
      graphics,
      oldTexture,
    );
    const { transform, texture, bbox, gridDim } = densityTextureData;

    const unitToCartn = Mat4.mul(Mat4(), transform, Mat4.fromScaling(Mat4(), gridDim));
    const cellDim = Mat4.getScaling(Vec3(), transform);
    const vol = DirectVolume.create(
      bbox,
      gridDim,
      transform,
      unitToCartn,
      cellDim,
      texture,
      stats,
      true,
      axisOrder,
      'byte',
      directVolume,
    );

    const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, densityTextureData.maxRadius);
    vol.setBoundingSphere(sphere);
    return vol;
  };

  return create(directVolume);
}

export const UnitsGaussianDensityVolumeParams = {
  ...UnitsDirectVolumeParams,
  jumpLength: JumpLengthParam,
  ...omitObjectKeys(GaussianDensityParams, ['floodfill']),
  ignoreHydrogens: PD.Boolean(false),
  ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
  includeParent: PD.Boolean(false, { isHidden: true }),
};
export type UnitsGaussianDensityVolumeParams = typeof UnitsGaussianDensityVolumeParams;

export function UnitsGaussianDensityVolumeVisual(materialId: number): UnitsVisual<UnitsGaussianDensityVolumeParams> {
  return UnitsDirectVolumeVisual<UnitsGaussianDensityVolumeParams>(
    {
      defaultProps: PD.getDefaultValues(UnitsGaussianDensityVolumeParams),
      createGeometry: createUnitsGaussianDensityVolume,
      // the density packs the group with the value
      signatureDefines: () => ({ dPackedGroup: [true] }),
      createLocationIterator: ElementIterator.fromGroup,
      getLoci: getElementLoci,
      eachLocation: eachElement,
      setUpdateState: (
        state: VisualUpdateState,
        newProps: PD.Values<GaussianDensityVolumeParams>,
        currentProps: PD.Values<GaussianDensityVolumeParams>,
      ) => {
        if (newProps.resolution !== currentProps.resolution) state.createGeometry = true;
        if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true;
        if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true;
        if (newProps.ignoreHydrogens !== currentProps.ignoreHydrogens) state.createGeometry = true;
        if (newProps.ignoreHydrogensVariant !== currentProps.ignoreHydrogensVariant) state.createGeometry = true;
        if (newProps.traceOnly !== currentProps.traceOnly) state.createGeometry = true;
        if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
      },
      dispose: (geometry: DirectVolume) => {
        geometry.gridTexture.ref.value.destroy();
      },
    },
    materialId,
  );
}
