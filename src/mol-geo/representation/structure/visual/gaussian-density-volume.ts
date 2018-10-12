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
import { GaussianDensityProps, GaussianDensityParams } from 'mol-model/structure/structure/unit/gaussian-density';
import { paramDefaultValues } from 'mol-view/parameter';
import { DirectVolume2d } from '../../../geometry/direct-volume/direct-volume';
import { ValueCell } from 'mol-util';
import { Vec3, Vec2 } from 'mol-math/linear-algebra';

async function createGaussianDensityVolume(ctx: RuntimeContext, unit: Unit, structure: Structure, props: GaussianDensityProps, directVolume?: DirectVolume2d): Promise<DirectVolume2d> {
    const p = { ...props, useGpu: true, ignoreCache: true }
    const { transform, renderTarget, bbox, gridDimension } = await unit.computeGaussianDensity(p, ctx)
    if (!renderTarget || !bbox || !gridDimension) throw new Error('missing renderTarget and/or boundingBox and/or gridDimension')

    if (directVolume) {
        ValueCell.update(directVolume.gridDimension, gridDimension)
        ValueCell.update(directVolume.gridTexture, renderTarget.image)
        ValueCell.update(directVolume.gridTextureDim, Vec2.set(directVolume.gridTextureDim.ref.value, renderTarget.width, renderTarget.height))
        ValueCell.update(directVolume.bboxMin, bbox.min)
        ValueCell.update(directVolume.bboxMax, bbox.max)
        ValueCell.update(directVolume.bboxSize, Vec3.sub(directVolume.bboxSize.ref.value, bbox.max, bbox.min))
        ValueCell.update(directVolume.transform, transform)
    } else {
        directVolume = {
            kind: 'direct-volume-2d' as 'direct-volume-2d',
            gridDimension: ValueCell.create(gridDimension),
            gridTexture: ValueCell.create(renderTarget.image),
            gridTextureDim: ValueCell.create(Vec2.create(renderTarget.width, renderTarget.height)),
            bboxMin: ValueCell.create(bbox.min),
            bboxMax: ValueCell.create(bbox.max),
            bboxSize: ValueCell.create(Vec3.sub(Vec3.zero(), bbox.max, bbox.min)),
            transform: ValueCell.create(transform),
        }
    }

    console.log('gridDimension', gridDimension)
    console.log('gridTextureDim', renderTarget.width, renderTarget.height)
    console.log('boundingBox', bbox)
    console.log('transform', transform)

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
            if (newProps.readSlices !== currentProps.readSlices) state.createGeometry = true
            if (newProps.ignoreCache !== currentProps.ignoreCache) state.createGeometry = true
        }
    })
}