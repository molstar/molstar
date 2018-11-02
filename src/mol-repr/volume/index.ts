/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RepresentationProps, Representation, Visual } from '..';
import { VolumeData } from 'mol-model/volume';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { paramDefaultValues } from 'mol-util/parameter';
import { Geometry } from 'mol-geo/geometry/geometry';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';

export interface VolumeVisual<P extends RepresentationProps = {}> extends Visual<VolumeData, P> { }

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export const VolumeParams = {
    ...Geometry.Params,
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
        label: 'Volume',
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