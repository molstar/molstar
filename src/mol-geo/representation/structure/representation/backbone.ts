/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation, UnitsRepresentation } from '..';
import { PickingId } from '../../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../../../util/marker-data';
import { PolymerBackboneVisual, DefaultPolymerBackboneProps } from '../visual/polymer-backbone-cylinder';

export const DefaultBackboneProps = {
    ...DefaultPolymerBackboneProps
}
export type BackboneProps = typeof DefaultBackboneProps

export type BackboneRepresentation = StructureRepresentation<BackboneProps>

export function BackboneRepresentation(): BackboneRepresentation {
    const traceRepr = UnitsRepresentation(PolymerBackboneVisual)

    let currentProps: BackboneProps
    return {
        get renderObjects() {
            return [ ...traceRepr.renderObjects ]
        },
        get props() {
            return { ...traceRepr.props }
        },
        createOrUpdate: (props: Partial<BackboneProps> = {}, structure?: Structure) => {
            currentProps = Object.assign({}, DefaultBackboneProps, props)
            return Task.create('BackboneRepresentation', async ctx => {
                await traceRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            return traceRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            traceRepr.mark(loci, action)
        },
        destroy() {
            traceRepr.destroy()
        }
    }
}