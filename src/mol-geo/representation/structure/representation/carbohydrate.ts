/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComplexRepresentation, StructureRepresentation } from '..';
import { PickingId } from '../../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci, isEmptyLoci } from 'mol-model/loci';
import { MarkerAction } from '../../../util/marker-data';
import { CarbohydrateSymbolVisual, DefaultCarbohydrateSymbolProps } from '../visual/carbohydrate-symbol-mesh';
import { CarbohydrateLinkVisual, DefaultCarbohydrateLinkProps } from '../visual/carbohydrate-link-cylinder';

export const DefaultCartoonProps = {
    ...DefaultCarbohydrateSymbolProps,
    ...DefaultCarbohydrateLinkProps
}
export type CarbohydrateProps = Partial<typeof DefaultCartoonProps>

export function CarbohydrateRepresentation(): StructureRepresentation<CarbohydrateProps> {
    const carbohydrateSymbolRepr = ComplexRepresentation(CarbohydrateSymbolVisual)
    const carbohydrateLinkRepr = ComplexRepresentation(CarbohydrateLinkVisual)

    return {
        get renderObjects() {
            return [ ...carbohydrateSymbolRepr.renderObjects, ...carbohydrateLinkRepr.renderObjects ]
        },
        get props() {
            return { ...carbohydrateSymbolRepr.props, ...carbohydrateLinkRepr.props }
        },
        create: (structure: Structure, props: CarbohydrateProps = {} as CarbohydrateProps) => {
            const p = Object.assign({}, DefaultCartoonProps, props)
            return Task.create('Creating CarbohydrateRepresentation', async ctx => {
                await carbohydrateSymbolRepr.create(structure, p).runInContext(ctx)
                await carbohydrateLinkRepr.create(structure, p).runInContext(ctx)
            })
        },
        update: (props: CarbohydrateProps) => {
            const p = Object.assign({}, props)
            return Task.create('Updating CarbohydrateRepresentation', async ctx => {
                await carbohydrateSymbolRepr.update(p).runInContext(ctx)
                await carbohydrateLinkRepr.update(p).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            const carbohydrateSymbolLoci = carbohydrateSymbolRepr.getLoci(pickingId)
            const carbohydrateLinkLoci = carbohydrateLinkRepr.getLoci(pickingId)
            return !isEmptyLoci(carbohydrateSymbolLoci) ? carbohydrateSymbolLoci
                : carbohydrateLinkLoci
        },
        mark: (loci: Loci, action: MarkerAction) => {
            carbohydrateSymbolRepr.mark(loci, action)
            carbohydrateLinkRepr.mark(loci, action)
        },
        destroy() {
            carbohydrateSymbolRepr.destroy()
            carbohydrateLinkRepr.destroy()
        }
    }
}