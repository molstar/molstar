/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { RepresentationProps, Representation, Visual } from '..';
import { VolumeData } from 'mol-model/volume';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';

export interface VolumeVisual<P extends RepresentationProps = {}> extends Visual<VolumeData, P> { }

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export function VolumeRepresentation<P>(visualCtor: (volumeData: VolumeData) => VolumeVisual<P>): VolumeRepresentation<P> {
    const renderObjects: RenderObject[] = []
    let _volumeData: VolumeData

    function create(volumeData: VolumeData, props: P = {} as P) {
        return Task.create('VolumeRepresentation.create', async ctx => {
            _volumeData = volumeData
            const visual = visualCtor(_volumeData)
            await visual.create(ctx, _volumeData, props)
            renderObjects.push(...visual.renderObjects)
        });
    }
    
    function update(props: P) {
        return Task.create('VolumeRepresentation.update', async ctx => {})
    }

    return {
        renderObjects,
        create,
        update,
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
        },
        destroy() {
            // TODO
        }
    }
}