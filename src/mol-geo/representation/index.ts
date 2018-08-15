/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task, RuntimeContext } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object'
import { PickingId } from '../util/picking';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../util/marker-data';

export interface RepresentationProps {}

export interface Representation<D, P extends RepresentationProps = {}> {
    readonly renderObjects: ReadonlyArray<RenderObject>
    readonly props: Readonly<P>
    create: (data: D, props?: Partial<P>) => Task<void>
    update: (props: Partial<P>) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
    destroy: () => void
}

export interface Visual<D, P extends RepresentationProps = {}> {
    readonly renderObject: RenderObject
    create: (ctx: RuntimeContext, data: D, props?: Partial<P>) => Promise<void>
    update: (ctx: RuntimeContext, props: Partial<P>) => Promise<boolean>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
    destroy: () => void
}