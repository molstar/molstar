/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation, UnitsRepresentation } from '..';
import { PickingId } from '../../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../../util/marker-data';
import { PolymerTraceVisual, DefaultPolymerTraceProps } from '../visual/polymer-trace-mesh';
import { PolymerGapVisual, DefaultPolymerGapProps } from '../visual/polymer-gap-cylinder';
import { NucleotideBlockVisual, DefaultNucleotideBlockProps } from '../visual/nucleotide-block-mesh';
import { PolymerDirectionVisual, DefaultPolymerDirectionProps } from '../visual/polymer-direction-wedge';

export const DefaultCartoonProps = {
    ...DefaultPolymerTraceProps,
    ...DefaultPolymerGapProps,
    ...DefaultNucleotideBlockProps,
    ...DefaultPolymerDirectionProps
}
export type CartoonProps = typeof DefaultCartoonProps

export function CartoonRepresentation(): StructureRepresentation<CartoonProps> {
    const traceRepr = UnitsRepresentation(PolymerTraceVisual)
    const gapRepr = UnitsRepresentation(PolymerGapVisual)
    const blockRepr = UnitsRepresentation(NucleotideBlockVisual)
    const directionRepr = UnitsRepresentation(PolymerDirectionVisual)

    let currentProps: CartoonProps
    return {
        get renderObjects() {
            return [ ...traceRepr.renderObjects, ...gapRepr.renderObjects,
                ...blockRepr.renderObjects // , ...directionRepr.renderObjects
            ]
        },
        get props() {
            return { ...traceRepr.props, ...gapRepr.props, ...blockRepr.props }
        },
        create: (structure: Structure, props: Partial<CartoonProps> = {}) => {
            currentProps = Object.assign({}, DefaultCartoonProps, props)
            return Task.create('Creating CartoonRepresentation', async ctx => {
                await traceRepr.create(structure, currentProps).runInContext(ctx)
                await gapRepr.create(structure, currentProps).runInContext(ctx)
                await blockRepr.create(structure, currentProps).runInContext(ctx)
                await directionRepr.create(structure, currentProps).runInContext(ctx)
            })
        },
        update: (props: Partial<CartoonProps>) => {
            currentProps = Object.assign(currentProps, props)
            return Task.create('Updating CartoonRepresentation', async ctx => {
                await traceRepr.update(currentProps).runInContext(ctx)
                await gapRepr.update(currentProps).runInContext(ctx)
                await blockRepr.update(currentProps).runInContext(ctx)
                await directionRepr.update(currentProps).runInContext(ctx)
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