/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Task } from 'mol-task'
import { RenderObject } from 'mol-gl/render-object'
import { PickingId } from '../util/picking';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../util/marker-data';

export interface RepresentationProps {}

export interface Representation<D, P extends RepresentationProps = {}> {
    renderObjects: ReadonlyArray<RenderObject>
    create: (data: D, props?: P) => Task<void>
    update: (props: P) => Task<void>
    getLoci: (pickingId: PickingId) => Loci
    mark: (loci: Loci, action: MarkerAction) => void
}