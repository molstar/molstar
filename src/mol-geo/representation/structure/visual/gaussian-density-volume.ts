/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Unit, Structure } from 'mol-model/structure';
import { UnitsVisual, VisualUpdateState } from '..';
import { RuntimeContext } from 'mol-task'
import { UnitsDirectVolumeVisual, UnitsDirectVolumeParams } from '../units-visual';
import { StructureElementIterator, getElementLoci, markElement } from './util/element';
import { GaussianDensityProps, GaussianDensityParams, computeUnitGaussianDensityTexture } from 'mol-model/structure/structure/unit/gaussian-density';
import { paramDefaultValues } from 'mol-view/parameter';
import { DirectVolume2d, DirectVolume3d } from '../../../geometry/direct-volume/direct-volume';

async function createGaussianDensityVolume(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianDensityProps, directVolume?: DirectVolume2d | DirectVolume3d): Promise<DirectVolume2d | DirectVolume3d> {
    const { webgl } = props
    if (webgl === undefined) throw new Error('createGaussianDensityVolume requires `webgl` in props')

    const p = { ...props, useGpu: true }
    const oldTexture = directVolume ? directVolume.gridTexture.ref.value : undefined
    const densityTextureData = await computeUnitGaussianDensityTexture(unit, p, oldTexture).runInContext(ctx)
    const { transform, texture, bbox, gridDimension } = densityTextureData

    directVolume = texture.depth === 0 ?
        DirectVolume2d.create(bbox, gridDimension, transform, texture, directVolume as DirectVolume2d) :
        DirectVolume3d.create(bbox, gridDimension, transform, texture, directVolume as DirectVolume3d)

    return directVolume;
}

export const GaussianDensityVolumeParams = {
    ...UnitsDirectVolumeParams,
    ...GaussianDensityParams,
}
export const DefaultGaussianDensityVolumeProps = paramDefaultValues(GaussianDensityVolumeParams)
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
            if (newProps.smoothness !== currentProps.smoothness) state.createGeometry = true
            if (newProps.useGpu !== currentProps.useGpu) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}