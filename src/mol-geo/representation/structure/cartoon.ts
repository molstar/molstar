/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation, StructureUnitsRepresentation } from '.';
import { PickingId } from '../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../util/marker-data';
import { PolymerTraceVisual, DefaultPolymerTraceProps } from './visual/polymer-trace-mesh';
import { PolymerGapVisual, DefaultPolymerGapProps } from './visual/polymer-gap-cylinder';
import { NucleotideBlockVisual, DefaultNucleotideBlockProps } from './visual/nucleotide-block-mesh';
import { PolymerDirectionVisual, DefaultPolymerDirectionProps } from './visual/polymer-direction-wedge';

export const DefaultCartoonProps = {
    ...DefaultPolymerTraceProps,
    ...DefaultPolymerGapProps,
    ...DefaultNucleotideBlockProps,
    ...DefaultPolymerDirectionProps
}
export type CartoonProps = Partial<typeof DefaultCartoonProps>

export function CartoonRepresentation(): StructureRepresentation<CartoonProps> {
    const traceRepr = StructureUnitsRepresentation(PolymerTraceVisual)
    const gapRepr = StructureUnitsRepresentation(PolymerGapVisual)
    const blockRepr = StructureUnitsRepresentation(NucleotideBlockVisual)
    const directionRepr = StructureUnitsRepresentation(PolymerDirectionVisual)

    return {
        get renderObjects() {
            return [ ...traceRepr.renderObjects, ...gapRepr.renderObjects, ...blockRepr.renderObjects, ...directionRepr.renderObjects ]
        },
        get props() {
            return { ...traceRepr.props, ...gapRepr.props, ...blockRepr.props, ...directionRepr.props }
        },
        create: (structure: Structure, props: CartoonProps = {} as CartoonProps) => {
            const p = Object.assign({}, DefaultCartoonProps, props)
            return Task.create('CartoonRepresentation', async ctx => {
                await traceRepr.create(structure, p).runInContext(ctx)
                await gapRepr.create(structure, p).runInContext(ctx)
                await blockRepr.create(structure, p).runInContext(ctx)
                await directionRepr.create(structure, p).runInContext(ctx)
            })
        },
        update: (props: CartoonProps) => {
            const p = Object.assign({}, props)
            return Task.create('Updating CartoonRepresentation', async ctx => {
                await traceRepr.update(p).runInContext(ctx)
                await gapRepr.update(p).runInContext(ctx)
                await blockRepr.update(p).runInContext(ctx)
                await directionRepr.update(p).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const traceLoci = traceRepr.getLoci(pickingId)
            const gapLoci = gapRepr.getLoci(pickingId)
            const blockLoci = blockRepr.getLoci(pickingId)
            const directionLoci = directionRepr.getLoci(pickingId)
            return !isEmptyLoci(traceLoci) ? traceLoci
                : !isEmptyLoci(gapLoci) ? gapLoci
                : !isEmptyLoci(blockLoci) ? blockLoci
                : directionLoci
        },
        mark: (loci: Loci, action: MarkerAction) => {
            traceRepr.mark(loci, action)
            gapRepr.mark(loci, action)
            blockRepr.mark(loci, action)
            directionRepr.mark(loci, action)
        },
        destroy() {
            traceRepr.destroy()
            gapRepr.destroy()
            blockRepr.destroy()
            directionRepr.destroy()
        }
    }
}