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
import { CrossLinkRestraintVisual, DefaultCrossLinkRestraintProps } from '../visual/cross-link-restraint-cylinder';
import { SizeThemeProps } from 'mol-view/theme/size';
import { getQualityProps } from '../../util';

export const DefaultDistanceRestraintProps = {
    ...DefaultCrossLinkRestraintProps,
    sizeTheme: { name: 'uniform', value: 0.25 } as SizeThemeProps,
}
export type DistanceRestraintProps = typeof DefaultDistanceRestraintProps

export type DistanceRestraintRepresentation = StructureRepresentation<DistanceRestraintProps>

export function DistanceRestraintRepresentation(): DistanceRestraintRepresentation {
    const crossLinkRepr = ComplexRepresentation('Cross-link restraint', CrossLinkRestraintVisual)

    let currentProps: DistanceRestraintProps
    return {
        label: 'Distance restraint',
        get renderObjects() {
            return [ ...crossLinkRepr.renderObjects ]
        },
        get props() {
            return { ...crossLinkRepr.props }
        },
        createOrUpdate: (props: Partial<DistanceRestraintProps> = {}, structure?: Structure) => {
            const qualityProps = getQualityProps(Object.assign({}, currentProps, props), structure)
            currentProps = Object.assign({}, DefaultDistanceRestraintProps, currentProps, props, qualityProps)
            return Task.create('DistanceRestraintRepresentation', async ctx => {
                await crossLinkRepr.createOrUpdate(currentProps, structure).runInContext(ctx)
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