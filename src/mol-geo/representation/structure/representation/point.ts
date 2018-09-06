/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { ElementPointVisual, DefaultElementPointProps } from '../visual/element-point';
import { StructureRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { MarkerAction } from '../../../util/marker-data';
import { Loci } from 'mol-model/loci';
import { PickingId } from '../../../util/picking';

export const DefaultPointProps = {
    ...DefaultElementPointProps,
}
export type PointProps = typeof DefaultPointProps

export type PointRepresentation = StructureRepresentation<PointProps>

export function PointRepresentation(): PointRepresentation {
    let currentProps: PointProps
    const pointRepr = UnitsRepresentation('Point', ElementPointVisual)
    return {
        label: 'Point',
        get renderObjects() {
            return [ ...pointRepr.renderObjects ]
        },
        get props() {
            return { ...pointRepr.props }
        },
        createOrUpdate: (props: Partial<PointProps> = {}, structure?: Structure) => {
            currentProps = Object.assign({}, DefaultPointProps, currentProps, props)
            return pointRepr.createOrUpdate(currentProps, structure)
        },
        getLoci: (pickingId: PickingId) => {
            return pointRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            pointRepr.mark(loci, action)
        },
        destroy() {
            pointRepr.destroy()
        }
    }
}