/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementPointVisual, ElementPointParams } from '../visual/element-point';
import { UnitsRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { Loci } from 'mol-model/loci';
import { paramDefaultValues } from 'mol-util/parameter';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { StructureRepresentation } from '../index';

export const PointParams = {
    ...ElementPointParams,
}
export const DefaultPointProps = paramDefaultValues(PointParams)
export type PointProps = typeof DefaultPointProps

export type PointRepresentation = StructureRepresentation<PointProps>

export function PointRepresentation(): PointRepresentation {
    let currentProps: PointProps
    const pointRepr = UnitsRepresentation('Point', ElementPointVisual)
    return {
        label: 'Point',
        params: PointParams,
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
            return pointRepr.mark(loci, action)
        },
        destroy() {
            pointRepr.destroy()
        }
    }
}