/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ComplexRepresentation, StructureRepresentation } from '..';
import { PickingId } from '../../../geometry/picking';
import { Structure } from 'mol-model/structure';
import { Task } from 'mol-task';
import { Loci } from 'mol-model/loci';
import { MarkerAction } from '../../../geometry/marker-data';
import { CrossLinkRestraintVisual, CrossLinkRestraintParams } from '../visual/cross-link-restraint-cylinder';
import { SizeThemeName, SizeThemeOptions } from 'mol-view/theme/size';
import { getQualityProps } from '../../util';
import { paramDefaultValues, SelectParam, NumberParam } from 'mol-view/parameter';

export const DistanceRestraintParams = {
    ...CrossLinkRestraintParams,
    sizeTheme: SelectParam<SizeThemeName>('Size Theme', '', 'uniform', SizeThemeOptions),
    sizeValue: NumberParam('Size Value', '', 0.25, 0, 0.05, 20),
}
export const DefaultDistanceRestraintProps = paramDefaultValues(DistanceRestraintParams)
export type DistanceRestraintProps = typeof DefaultDistanceRestraintProps

export type DistanceRestraintRepresentation = StructureRepresentation<DistanceRestraintProps>

export function DistanceRestraintRepresentation(): DistanceRestraintRepresentation {
    const crossLinkRepr = ComplexRepresentation('Cross-link restraint', CrossLinkRestraintVisual)

    let currentProps: DistanceRestraintProps
    return {
        label: 'Distance restraint',
        params: DistanceRestraintParams,
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
            return crossLinkRepr.mark(loci, action)
        },
        destroy() {
            crossLinkRepr.destroy()
        }
    }
}