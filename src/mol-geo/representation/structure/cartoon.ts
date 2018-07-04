/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation, StructureUnitsRepresentation } from '.';
import { PickingId } from '../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';
import { PolymerTraceVisual, DefaultPolymerTraceProps } from './visual/polymer-trace-mesh';

export const DefaultCartoonProps = {
    ...DefaultPolymerTraceProps
}
export type CartoonProps = Partial<typeof DefaultCartoonProps>

export function CartoonRepresentation(): StructureRepresentation<CartoonProps> {
    const traceRepr = StructureUnitsRepresentation(PolymerTraceVisual)

    return {
        get renderObjects() {
            return [ ...traceRepr.renderObjects ]
        },
        get props() {
            return { ...traceRepr.props }
        },
        create: (structure: Structure, props: CartoonProps = {} as CartoonProps) => {
            const p = Object.assign({}, DefaultCartoonProps, props)
            return Task.create('CartoonRepresentation', async ctx => {
                await traceRepr.create(structure, p).runInContext(ctx)
            })
        },
        update: (props: CartoonProps) => {
            const p = Object.assign({}, props)
            return Task.create('Updating CartoonRepresentation', async ctx => {
                await traceRepr.update(p).runInContext(ctx)
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