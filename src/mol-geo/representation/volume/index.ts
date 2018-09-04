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
import { DefaultBaseProps } from '../util';

export interface VolumeVisual<P extends RepresentationProps = {}> extends Visual<VolumeData, P> { }

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export const DefaultVolumeProps = {
    ...DefaultBaseProps
}
export type VolumeProps = typeof DefaultVolumeProps

export function VolumeRepresentation<P extends VolumeProps>(visualCtor: (volumeData: VolumeData) => VolumeVisual<P>): VolumeRepresentation<P> {
    const renderObjects: RenderObject[] = []
    let _volumeData: VolumeData
    let _props: P

    function createOrUpdate(props: Partial<P> = {}, volumeData?: VolumeData) {
        _props = Object.assign({}, DefaultVolumeProps, _props, props)
        return Task.create('VolumeRepresentation.create', async ctx => {
            if (volumeData) {
                _volumeData = volumeData
                const visual = visualCtor(_volumeData)
                await visual.createOrUpdate(ctx, props, _volumeData)
                if (visual.renderObject) renderObjects.push(visual.renderObject)
            } else {
                throw new Error('missing volumeData')
            }
        });
    }

    return {
        label: 'Volume mesh',
        get renderObjects () { return renderObjects },
        get props () { return _props },
        createOrUpdate,
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