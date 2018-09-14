/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StructureRepresentation, UnitsRepresentation } from '..';
import { PickingId } from '../../../geometry/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../../geometry/marker-data';
import { PolymerTraceVisual, DefaultPolymerTraceProps } from '../visual/polymer-trace-mesh';
import { PolymerGapVisual, DefaultPolymerGapProps } from '../visual/polymer-gap-cylinder';
import { NucleotideBlockVisual, DefaultNucleotideBlockProps } from '../visual/nucleotide-block-mesh';
import { SizeThemeProps } from 'mol-view/theme/size';
import { getQualityProps } from '../../util';
// import { PolymerDirectionVisual, DefaultPolymerDirectionProps } from '../visual/polymer-direction-wedge';

export const DefaultCartoonProps = {
    ...DefaultPolymerTraceProps,
    ...DefaultPolymerGapProps,
    ...DefaultNucleotideBlockProps,
    // ...DefaultPolymerDirectionProps,

    sizeTheme: { name: 'uniform', value: 0.2 } as SizeThemeProps,
}
export type CartoonProps = typeof DefaultCartoonProps

export type CartoonRepresentation = StructureRepresentation<CartoonProps>

export function CartoonRepresentation(): CartoonRepresentation {
    const traceRepr = UnitsRepresentation('Polymer trace mesh', PolymerTraceVisual)
    const gapRepr = UnitsRepresentation('Polymer gap cylinder', PolymerGapVisual)
    const blockRepr = UnitsRepresentation('Nucleotide block mesh', NucleotideBlockVisual)
    // const directionRepr = UnitsRepresentation('Polymer direction wedge', PolymerDirectionVisual)

    let currentProps: CartoonProps
    return {
        label: 'Cartoon',
        get renderObjects() {
            return [ ...traceRepr.renderObjects, ...gapRepr.renderObjects,
                ...blockRepr.renderObjects // , ...directionRepr.renderObjects
            ]
        },
        get props() {
            return { ...traceRepr.props, ...gapRepr.props, ...blockRepr.props }
        },
        createOrUpdate: (props: Partial<CartoonProps> = {}, structure?: Structure) => {
            const qualityProps = getQualityProps(Object.assign({}, currentProps, props), structure)
            currentProps = Object.assign({}, DefaultCartoonProps, currentProps, props, qualityProps)
            return Task.create('Creating CartoonRepresentation', async ctx => {
                await traceRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
                await gapRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
                await blockRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
                // await directionRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const traceLoci = traceRepr.getLoci(pickingId)
            const gapLoci = gapRepr.getLoci(pickingId)
            const blockLoci = blockRepr.getLoci(pickingId)
            // const directionLoci = directionRepr.getLoci(pickingId)
            return !isEmptyLoci(traceLoci) ? traceLoci
                : !isEmptyLoci(gapLoci) ? gapLoci
                : blockLoci
                // : !isEmptyLoci(blockLoci) ? blockLoci
                // : directionLoci
        },
        mark: (loci: Loci, action: MarkerAction) => {
            const markTrace = traceRepr.mark(loci, action)
            const markGap = gapRepr.mark(loci, action)
            const markBlock = blockRepr.mark(loci, action)
            // const markDirection = directionRepr.mark(loci, action)
            return markTrace || markGap || markBlock // \\ markDirection
        },
        destroy() {
            traceRepr.destroy()
            gapRepr.destroy()
            blockRepr.destroy()
            // directionRepr.destroy()
        }
    }
}