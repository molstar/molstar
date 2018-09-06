/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { UnitsRepresentation } from '..';
import { ElementSphereVisual, DefaultElementSphereProps } from '../visual/element-sphere';
import { StructureRepresentation } from '../units-representation';
import { Structure } from 'mol-model/structure';
import { PickingId } from '../../../util/picking';
import { MarkerAction } from '../../../util/marker-data';
import { Loci } from 'mol-model/loci';
import { getQualityProps } from '../../util';

export const DefaultSpacefillProps = {
    ...DefaultElementSphereProps
}
export type SpacefillProps = typeof DefaultSpacefillProps

export type SpacefillRepresentation = StructureRepresentation<SpacefillProps>

export function SpacefillRepresentation(): SpacefillRepresentation {
    let currentProps: SpacefillProps
    const sphereRepr = UnitsRepresentation('Sphere mesh', ElementSphereVisual)
    return {
        label: 'Spacefill',
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
            sphereRepr.mark(loci, action)
        },
        destroy() {
            sphereRepr.destroy()
        }
    }
}