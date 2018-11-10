/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual } from '../representation';
import { VisualUpdateState } from '../../util';
import { UnitsDirectVolumeVisual, UnitsDirectVolumeParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { GaussianDensityProps, GaussianDensityParams, computeUnitGaussianDensityTexture } from 'mol-model/structure/structure/unit/gaussian-density';
import { ParamDefinition as PD } from 'mol-util/param-definition';
import { DirectVolume } from 'mol-geo/geometry/direct-volume/direct-volume';
import { VisualContext } from 'mol-repr/representation';
import { Theme } from 'mol-theme/theme';

async function createGaussianDensityVolume(ctx: VisualContext, unit: Unit, structure: Structure, theme: Theme, props: GaussianDensityProps, directVolume?: DirectVolume): Promise<DirectVolume> {
    const { runtime, webgl } = ctx
    if (webgl === undefined) throw new Error('createGaussianDensityVolume requires `webgl` object in VisualContext')

    const p = { ...props, useGpu: true }
    const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined
    const densityTextureData = await computeUnitGaussianDensityTexture(unit, p, webgl, oldTexture).runInContext(runtime)
    const { transform, texture, bbox, gridDimension } = densityTextureData

    return DirectVolume.create(bbox, gridDimension, transform, texture, directVolume)
}

export const GaussianDensityVolumeParams = {
    ...UnitsDirectVolumeParams,
    ...GaussianDensityParams,
}
export const DefaultGaussianDensityVolumeProps = PD.getDefaultValues(GaussianDensityVolumeParams)
export type GaussianDensityVolumeProps = typeof DefaultGaussianDensityVolumeProps

export function GaussianDensityVolumeVisual(): UnitsVisual<GaussianDensityVolumeProps> {
    return UnitsDirectVolumeVisual<GaussianDensityVolumeProps>({
        defaultProps: DefaultGaussianDensityVolumeProps,
        createGeometry: createGaussianDensityVolume,
        createLocationIterator: StructureElementIterator.fromGroup,
        getLoci: getElementLoci,
        mark: markElement,
        setUpdateState: (state: VisualUpdateState, newProps: GaussianDensityVolumeProps, currentProps: GaussianDensityVolumeProps) => {
            if (newProps.resolution !== currentProps.resolution) state.createGeometry = true
            if (newProps.radiusOffset !== currentProps.radiusOffset) state.createGeometry = true
            if (newProps.smoothness !== currentProps.smoothness) {
                state.createGeometry = true
                newProps.isoValueAbsolute = Math.exp(-newProps.smoothness)
            }
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}