/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RepresentationProps, Representation, Visual } from '..';
import { VolumeData } from 'mol-model/volume';
import { PickingId } from '../../geometry/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../geometry/marker-data';
import { Geometry } from '../../geometry/geometry';
import { paramDefaultValues } from 'mol-view/parameter';
import { IsosurfaceParams } from './isosurface';

export interface VolumeVisual<P extends RepresentationProps = {}> extends Visual<VolumeData, P> { }

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export const VolumeParams = {
    ...Geometry.Params,
    ...IsosurfaceParams
}
export const DefaultVolumeProps = paramDefaultValues(VolumeParams)
export type VolumeProps = typeof DefaultVolumeProps

export function VolumeRepresentation<P extends VolumeProps>(visualCtor: (volumeData: VolumeData) => VolumeVisual<P>): VolumeRepresentation<P> {
    let visual: VolumeVisual<any>
    let _props: P
    let busy = false

    function createOrUpdate(props: Partial<P> = {}, volumeData?: VolumeData) {
        _props = Object.assign({}, DefaultVolumeProps, _props, props)
        return Task.create('VolumeRepresentation.create', async ctx => {
            // TODO queue it somehow
            if (busy) return

            if (!visual && !volumeData) {
                throw new Error('volumeData missing')
            } else if (volumeData && !visual) {
                busy = true
                visual = visualCtor(volumeData)
                await visual.createOrUpdate(ctx, props, volumeData)
                busy = false
            } else {
                busy = true
                await visual.createOrUpdate(ctx, props, volumeData)
                busy = false
            }
        });
    }

    return {
        label: 'Volume mesh',
        params: VolumeParams,
        get renderObjects() {
            return visual && visual.renderObject ? [ visual.renderObject ] : []
        },
        get props () { return _props },
        createOrUpdate,
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
            return false
        },
        destroy() {
            // TODO
        }
    }
}