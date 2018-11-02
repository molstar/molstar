/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ElementSphereVisual, ElementSphereParams } from '../visual/element-sphere';
import { UnitsRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { Loci } from 'mol-model/loci';
import { getQualityProps } from '../../util';
import { paramDefaultValues } from 'mol-util/parameter';
import { PickingId } from 'mol-geo/geometry/picking';
import { MarkerAction } from 'mol-geo/geometry/marker-data';
import { StructureRepresentation } from '../index';

export const SpacefillParams = {
    ...ElementSphereParams
}
export const DefaultSpacefillProps = paramDefaultValues(SpacefillParams)
export type SpacefillProps = typeof DefaultSpacefillProps

export type SpacefillRepresentation = StructureRepresentation<SpacefillProps>

export function SpacefillRepresentation(): SpacefillRepresentation {
    let currentProps: SpacefillProps
    const sphereRepr = UnitsRepresentation('Sphere mesh', ElementSphereVisual)
    return {
        label: 'Spacefill',
        params: SpacefillParams,
        get renderObjects() {
            return [ ...sphereRepr.renderObjects ]
        },
        get props() {
            return { ...sphereRepr.props }
        },
        createOrUpdate: (props: Partial<SpacefillProps> = {}, structure?: Structure) => {
            const qualityProps = getQualityProps(Object.assign({}, currentProps, props), structure)
            currentProps = Object.assign({}, DefaultSpacefillProps, currentProps, props, qualityProps)
            return sphereRepr.createOrUpdate(currentProps, structure)
        },
        getLoci: (pickingId: PickingId) => {
            return sphereRepr.getLoci(pickingId)
        },
        mark: (loci: Loci, action: MarkerAction) => {
            return sphereRepr.mark(loci, action)
        },
        destroy() {
            sphereRepr.destroy()
        }
    }
}