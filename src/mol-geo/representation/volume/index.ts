/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/scene';
import { RepresentationProps, Representation } from '..';
import { VolumeData } from 'mol-model/volume';

export interface VolumeElementRepresentation<P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (volumeData: VolumeData, props: P) => Task<void>
    update: (props: P) => Task<boolean>
}

export interface VolumeRepresentation<P extends RepresentationProps = {}> extends Representation<VolumeData, P> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (volumeData: VolumeData, props?: P) => Task<void>
    update: (props: P) => Task<void>
}

export function VolumeRepresentation<P>(reprCtor: () => VolumeElementRepresentation<P>): VolumeRepresentation<P> {
    const renderObjects: RenderObject[] = []

    return {
        renderObjects,
        create(volumeData: VolumeData, props: P = {} as P) {
            return Task.create('VolumeRepresentation.create', async ctx => {
                const repr = reprCtor()
                await ctx.update({ message: 'Building volume representation...', current: 0, max: 1 });
                await ctx.runChild(repr.create(volumeData, props));
                renderObjects.push(...repr.renderObjects)
            });
        },
        update(props: P) {
            return Task.create('VolumeRepresentation.update', async ctx => {})
        }
    }
}