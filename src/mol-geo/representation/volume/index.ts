/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object';
import { RepresentationProps, Representation } from '..';
import { VolumeData } from 'mol-model/volume';
import { PickingId } from '../../util/picking';
import { Loci, EmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';

export interface VolumeElementRepresentation<P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (volumeData: VolumeData, props: P) => Task<void>
    update: (props: P) => Task<boolean>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
}

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> { }

export function VolumeRepresentation<P>(reprCtor: () => VolumeElementRepresentation<P>): VolumeRepresentation<P> {
    const renderObjects: RenderObject[] = []

    return {
        renderObjects,
        create(volumeData: VolumeData, props: P = {} as P) {
            return Task.create('VolumeRepresentation.create', async ctx => {
                const repr = reprCtor()
                await repr.create(volumeData, props).runAsChild(ctx, { message: 'Building volume representation...', current: 0, max: 1 });
                renderObjects.push(...repr.renderObjects)
            });
        },
        update(props: P) {
            return Task.create('VolumeRepresentation.update', async ctx => {})
        },
        getLoci(pickingId: PickingId) {
            // TODO
            return EmptyLoci
        },
        mark(loci: Loci, action: MarkerAction) {
            // TODO
        }
    }
}