/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Structure } from 'mol-model/structure';
import { VisualUpdateState } from '../../util';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { VisualContext } from 'mol-repr/visual';
import { Theme } from 'mol-theme/theme';
import { ComplexVisual, ComplexDirectVolumeVisual, ComplexDirectVolumeParams } from '../complex-visual';
import { EmptyLoci } from 'mol-model/loci';
import { NullLocation } from 'mol-model/location';
import { LocationIterator } from 'mol-geo/util/location-iterator';
import { computeStructureGaussianDensityTexture, GaussianDensityTextureProps, GaussianDensityTextureParams } from './util/gaussian';

async function createGaussianDensityVolume(ctx: VisualContext, structure: Structure, theme: Theme, props: GaussianDensityTextureProps, directVolume?: DirectVolume): Promise<DirectVolume> {
    const { runtime, webgl } = ctx
    if (webgl === undefined) throw new Error('createGaussianDensityVolume requires `webgl` object in VisualContext')

    const p = { ...props, useGpu: true }
    const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined
    const densityTextureData = await computeStructureGaussianDensityTexture(structure, p, webgl, oldTexture).runInContext(runtime)
    const { transform, texture, bbox, gridDimension } = densityTextureData

    return DirectVolume.create(bbox, gridDimension, transform, texture, directVolume)
}

export const GaussianDensityVolumeParams = {
    ...ComplexDirectVolumeParams,
    ...GaussianDensityTextureParams
}
export type GaussianDensityVolumeParams = typeof GaussianDensityVolumeParams

export function GaussianDensityVolumeVisual(): ComplexVisual<GaussianDensityVolumeParams> {
    return ComplexDirectVolumeVisual<GaussianDensityVolumeParams>({
        defaultProps: PD.getDefaultValues(GaussianDensityVolumeParams),
        createGeometry: createGaussianDensityVolume,
        createLocationIterator: (structure: Structure) => LocationIterator(structure.elementCount, 1, () => NullLocation),
        getLoci: () => EmptyLoci, // TODO
        eachLocation: () => false, // TODO
        setUpdateState: (state: VisualUpdateState, newProps: PD.Values<GaussianDensityVolumeParams>, currentProps: PD.Values<GaussianDensityVolumeParams>) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) {
                state.createGeometry = true
                newProps.isoValueNorm = Math.exp(-newProps.smoothness)
            }
        }
    })
}