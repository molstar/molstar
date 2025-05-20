/**
 * Copyright (c) 2018-2025 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure, Unit } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { GaussianDensityProps, computeStructureGaussianDensityTexture, computeUnitGaussianDensityTexture, GaussianDensityParams } from './util/gaussian';
import { DirectVolume } from '../../../mol-geo/geometry/direct-volume/direct-volume';
import { ComplexDirectVolumeParams, ComplexVisual, ComplexDirectVolumeVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { eachElement, eachSerialElement, ElementIterator, getElementLoci, getSerialElementLoci } from './util/element';
import { Sphere3D } from '../../../mol-math/geometry';
import { UnitsDirectVolumeParams, UnitsVisual, UnitsDirectVolumeVisual } from '../units-visual';

function createGaussianDensityVolume(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityProps, directVolume?: DirectVolume): DirectVolume {
    const { webgl } = ctx;
    if (!webgl) {
        // gpu gaussian density also needs blendMinMax but there is no fallback here so
        // we allow it here with the results that there is no group id assignment and
        // hence no group-based coloring or picking
        throw new Error('GaussianDensityVolume requires `webgl`');
    }

    const axisOrder = Vec3.create(0, 1, 2);
    const stats = { min: 0, max: 1, mean: 0.04, sigma: 0.01 };

    const create = (directVolume?: DirectVolume) => {
        const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined;
        const densityTextureData = computeStructureGaussianDensityTexture(structure, theme.size, props, webgl, oldTexture);
        const { transform, texture, bbox, gridDim } = densityTextureData;

        const unitToCartn = Mat4.mul(Mat4(), transform, Mat4.fromScaling(Mat4(), gridDim));
        const cellDim = Mat4.getScaling(Vec3(), transform);

        const vol = DirectVolume.create(bbox, gridDim, transform, unitToCartn, cellDim, texture, stats, true, axisOrder, 'byte', directVolume);

        const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, densityTextureData.maxRadius);
        vol.setBoundingSphere(sphere);
        return vol;
    };

    const vol = create(directVolume);
    vol.meta.reset = () => {
        create(vol);
    };

    return vol;
}

export const GaussianDensityVolumeParams = {
    ...ComplexDirectVolumeParams,
    ...GaussianDensityParams,
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    includeParent: PD.Boolean(false, { isHidden: true }),
};
export type GaussianDensityVolumeParams = typeof GaussianDensityVolumeParams

export function GaussianDensityVolumeVisual(materialId: number): ComplexVisual<GaussianDensityVolumeParams> {
    return ComplexDirectVolumeVisual<GaussianDensityVolumeParams>({
        defaultProps: PD.getDefaultValues(GaussianDensityVolumeParams),
        createGeometry: createGaussianDensityVolume,
        createLocationIterator: ElementIterator.fromStructure,
        getLoci: getSerialElementLoci,
        eachLocation: eachSerialElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianDensityVolumeParams>, currentProps: PD.Values<GaussianDensityVolumeParams>) => {
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
        }
    }, materialId);
}

//

function createUnitsGaussianDensityVolume(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, directVolume?: DirectVolume): DirectVolume {
    const { webgl } = ctx;
    if (!webgl) {
        // gpu gaussian density also needs blendMinMax but there is no fallback here so
        // we allow it here with the results that there is no group id assignment and
        // hence no group-based coloring or picking
        throw new Error('GaussianDensityVolume requires `webgl`');
    }

    const axisOrder = Vec3.create(0, 1, 2);
    const stats = { min: 0, max: 1, mean: 0.04, sigma: 0.01 };

    const create = (directVolume?: DirectVolume) => {
        const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined;
        const densityTextureData = computeUnitGaussianDensityTexture(structure, unit, theme.size, props, webgl, oldTexture);
        const { transform, texture, bbox, gridDim } = densityTextureData;

        const unitToCartn = Mat4.mul(Mat4(), transform, Mat4.fromScaling(Mat4(), gridDim));
        const cellDim = Mat4.getScaling(Vec3(), transform);
        const vol = DirectVolume.create(bbox, gridDim, transform, unitToCartn, cellDim, texture, stats, true, axisOrder, 'byte', directVolume);

        const sphere = Sphere3D.expand(Sphere3D(), unit.boundary.sphere, densityTextureData.maxRadius);
        vol.setBoundingSphere(sphere);
        return vol;
    };

    const vol = create(directVolume);
    vol.meta.reset = () => {
        create(vol);
    };

    return vol;
}

export const UnitsGaussianDensityVolumeParams = {
    ...UnitsDirectVolumeParams,
    ...GaussianDensityParams,
    ignoreHydrogens: PD.Boolean(false),
    ignoreHydrogensVariant: PD.Select('all', PD.arrayToOptions(['all', 'non-polar'] as const)),
    includeParent: PD.Boolean(false, { isHidden: true }),
};
export type UnitsGaussianDensityVolumeParams = typeof UnitsGaussianDensityVolumeParams

export function UnitsGaussianDensityVolumeVisual(materialId: number): UnitsVisual<UnitsGaussianDensityVolumeParams> {
    return UnitsDirectVolumeVisual<UnitsGaussianDensityVolumeParams>({
        defaultProps: PD.getDefaultValues(UnitsGaussianDensityVolumeParams),
        createGeometry: createUnitsGaussianDensityVolume,
        createLocationIterator: ElementIterator.fromGroup,
        getLoci: getElementLoci,
        eachLocation: eachElement,
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianDensityVolumeParams>, currentProps: PD.Values<GaussianDensityVolumeParams>) => {
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
        }
    }, materialId);
}