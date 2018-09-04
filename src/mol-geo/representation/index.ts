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
    readonly label: string
    readonly renderObjects: ReadonlyArray<RenderObject>
    readonly props: Readonly<P>
    createOrUpdate: (props?: Partial<P>, data?: D) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
    destroy: () => void
}

export interface Visual<D, P extends RepresentationProps = {}> {
    readonly renderObject: RenderObject | undefined
    createOrUpdate: (ctx: RuntimeContext, props?: Partial<P>, data?: D) => Promise<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
    destroy: () => void
}