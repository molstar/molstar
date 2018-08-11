/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComplexRepresentation, StructureRepresentation } from '..';
import { PickingId } from '../../../util/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../../../util/marker-data';
import { SizeTheme } from '../../../theme';
import { CrossLinkRestraintVisual, DefaultCrossLinkRestraintProps } from '../visual/cross-link-restraint-cylinder';

export const DefaultDistanceRestraintProps = {
    ...DefaultCrossLinkRestraintProps,
    sizeTheme: { name: 'uniform', value: 0.25 } as SizeTheme,
}
export type DistanceRestraintProps = typeof DefaultDistanceRestraintProps

export function DistanceRestraintRepresentation(): StructureRepresentation<DistanceRestraintProps> {
    const crossLinkRepr = ComplexRepresentation(CrossLinkRestraintVisual)

    let currentProps: DistanceRestraintProps
    return {
        get renderObjects() {
            return [ ...crossLinkRepr.renderObjects ]
        },
        get props() {
            return { ...crossLinkRepr.props }
        },
        create: (structure: Structure, props: Partial<DistanceRestraintProps> = {}) => {
            currentProps = Object.assign({}, DefaultDistanceRestraintProps, props)
            return Task.create('DistanceRestraintRepresentation', async ctx => {
                await crossLinkRepr.create(structure, currentProps).runInContext(ctx)
            })
        },
        update: (props: Partial<DistanceRestraintProps>) => {
            currentProps = Object.assign(currentProps, props)
            return Task.create('Updating DistanceRestraintRepresentation', async ctx => {
                await crossLinkRepr.update(currentProps).runInContext(ctx)
            })
        },
        getLoci: (pickingId: PickingId) => {
            return crossLinkRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            crossLinkRepr.mark(loci, action)
        },
        destroy() {
            crossLinkRepr.destroy()
        }
    }
}