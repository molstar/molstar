/**
 * Copyright (c) 2018-2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ParamDefinition as PD } from '../../../mol-util/param-definition';
import { VisualContext } from '../../visual';
import { Structure } from '../../../mol-model/structure';
import { Theme } from '../../../mol-theme/theme';
import { GaussianDensityTextureProps, computeStructureGaussianDensityTexture, GaussianDensityTextureParams } from './util/gaussian';
import { DirectVolume } from '../../../mol-geo/geometry/direct-volume/direct-volume';
import { ComplexDirectVolumeParams, ComplexVisual, ComplexDirectVolumeVisual } from '../complex-visual';
import { VisualUpdateState } from '../../util';
import { Mat4, Vec3 } from '../../../mol-math/linear-algebra';
import { eachSerialElement, ElementIterator, getSerialElementLoci } from './util/element';
import { Sphere3D } from '../../../mol-math/geometry';

async function createGaussianDensityVolume(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityTextureProps, directVolume?: DirectVolume): Promise<DirectVolume> {
    const { runtime, webgl } = ctx;
    if (webgl === undefined) throw new Error('createGaussianDensityVolume requires `webgl` object in VisualContext');

    const p = { ...props, useGpu: true };
    const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined;
    const densityTextureData = await computeStructureGaussianDensityTexture(structure, p, webgl, oldTexture).runInContext(runtime);
    const { transform, texture, bbox, gridDim } = densityTextureData;
    const stats = { min: 0, max: 1, mean: 0.04, sigma: 0.01 };

    const unitToCartn = Mat4.mul(Mat4(), transform, Mat4.fromScaling(Mat4(), gridDim));
    const cellDim = Vec3.create(1, 1, 1);

    const vol = DirectVolume.create(bbox, gridDim, transform, unitToCartn, cellDim, texture, stats, true, directVolume);

    const sphere = Sphere3D.expand(Sphere3D(), structure.boundary.sphere, props.radiusOffset);
    vol.setBoundingSphere(sphere);

    return vol;
}

export const GaussianDensityVolumeParams = {
    ...ComplexDirectVolumeParams,
    ...GaussianDensityTextureParams,
    ignoreHydrogens: PD.Boolean(false),
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
            if (newProps.includeParent !== currentProps.includeParent) state.createGeometry = true;
        }
    }, materialId);
}